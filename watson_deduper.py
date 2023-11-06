#!/usr/bin/env python
#Author: André Watson, November 2023

import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description = "Deduplicate aligned single-end reads from a SAM file. Assumes all duplicates are sequential: \
                                     be sure to run \"watson_deduper_setup.py\" followed by sorting the intermediate file produced.")
    parser.add_argument("-f", "--file", help="Absolute path to adjusted, sorted SAM file", required = True)
    parser.add_argument("-o", "--outfile", help="Absolute path to output deduplicated SAM file", required = True)
    parser.add_argument("-u", "--umi", help = "File containing list of valid UMIs", required = True)
    parser.add_argument("-m", "--mode", help = "Which duplicate to keep. Options include \"first\", \"last\", and \"best\" (highest MAPQ)", default="first")
    return parser.parse_args()


def get_umi(alignment: list[str]) -> str:
    '''
    Extracts and returns an UMI from the last colon-separated element of the QNAME field of an alignment line.
    '''
    return alignment[0].split(":")[-1]


def restore_pos(alignment: list[str]):
    '''
    The reciprocal operation of adjust_pos(). Restores the original `POS` field of an alignment line.
    '''
    old_pos = alignment.pop().split(":")[-1]
    alignment[3] = old_pos

#main script

#read in arguments
args = get_args()
file : str = args.file
outfile: str = args.outfile
umi_file: str = args.umi
mode: str = args.mode

#define valid UMIs
umi_set: set[str] = set()
with open(umi_file,"rt") as umifile:
    for line in umifile:
        umi_set.add(line.strip())

#set up duplicate tracker and convenience variables
comparing_mol: list[str] = []
bad_umi_ctr: int = 0
bad_mapping_ctr: int = 0
duplicate_ctr: int = 0
non_dupe_ctr: int = 0
first_alignment: bool = True
comparing_rev: bool = False
comparing_umi: str = ""

#open files
out = open(outfile, "wt")
un_deduped = open(file, "rt")

#loop over input file
while(True):
    line = un_deduped.readline().strip()
    #end of file:
    if line == "":
        #write out our last stored non-header line
        if len(comparing_mol) != 0:
            out.write("\t".join(comparing_mol) + "\n")
        break
    
    #header lines:
    if line.startswith("@"):
        out.write(line + "\n")
        continue
    
    #alignment lines:
    #split alignment line into fields
    alignment = line.split("\t")
    #check mapping
    if int(alignment[1]) & 4 == 4:
        #unmapped
        bad_mapping_ctr += 1
        continue
    
    #check UMI validity
    umi: str = get_umi(alignment)
    if umi not in umi_set:
        bad_umi_ctr += 1
        continue

    #grab first alignment line and store it for future comparisons
    if first_alignment:
        non_dupe_ctr += 1
        first_alignment = False
        comparing_mol = alignment
        comparing_rev = int(comparing_mol[1]) & 16 == 16 
        comparing_umi = umi
        continue

    
    #get strandedness of most recently read alignment line
    rev: bool = int(alignment[1]) & 16 == 16

    if(rev == comparing_rev
       and umi == comparing_umi
       and alignment[2] == comparing_mol[2]
       and alignment[3] == comparing_mol[3]):
        #we have a duplicate
        duplicate_ctr += 1
        match mode:
            case "first":
                #ignore all duplicates after the first, don't write them out
                pass
            case "last":
                #always replace stored alignment line with the most recent duplicate
                comparing_mol = alignment
                #strandedness and umi convenience variables don't need to be updated, since they're guaranteed to remain the same
            case "best":
                #replace stored alignment line if the current duplicate has a known, superior mapping quality 
                if int(alignment[4] != 255) and int(alignment[4]) > int(comparing_mol[4]):
                    comparing_mol = alignment
            case _:
                raise ValueError("Invalid option given for \"mode\"")
    else:
        #non-duplicate
        non_dupe_ctr += 1
        #write out stored alignment line
        restore_pos(comparing_mol)
        out.write("\t".join(comparing_mol) + "\n")
        #store new non-duplicate alignment line to check /it/ for duplicates
        comparing_mol = alignment
        comparing_rev = int(comparing_mol[1]) & 16 == 16
        comparing_umi = umi

#close files
out.close()
un_deduped.close()

print(f"Finished deduplicating file.\
        \nRemoved {bad_mapping_ctr} reads due to their being unmapped.\
        \nRemoved {bad_umi_ctr} reads due to invalid UMIs.\
        \nRemoved {duplicate_ctr} duplicate reads.\
        \nReads remaining: {non_dupe_ctr} unique reads.")
    
