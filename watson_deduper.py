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
same_pos_alignments: dict[tuple[str,bool],list[str]] = {}

#open files
out = open(outfile, "wt")
un_deduped = open(file, "rt")

#loop over input file
while(True):
    line = un_deduped.readline().strip()
    #end of file:
    if line == "":
        #write out our last stored non-header line
        for al in same_pos_alignments.values():
            non_dupe_ctr += 1
            restore_pos(al)
            out.write("\t".join(al) + "\n")
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

    #get strandedness of most recently read alignment line
    rev: bool = int(alignment[1]) & 16 == 16

    #grab first alignment line and store it for future comparisons
    if first_alignment:
        first_alignment = False
        comparing_mol = alignment
        same_pos_alignments[(umi, rev)] = comparing_mol
        continue


    if(alignment[2] == comparing_mol[2]
       and alignment[3] == comparing_mol[3]):
        #we are in a stretch of potential duplicates (same position and chromosome)
        if (umi, rev) not in same_pos_alignments:
            same_pos_alignments[(umi,rev)] = alignment
        else:
            duplicate_ctr += 1
            match mode:
                case "first":
                    pass
                case "last":
                    same_pos_alignments[(umi,rev)] = alignment
                case "best":
                    if int(alignment[4] != 255 and int(alignment[4]) > int(same_pos_alignments[(umi,rev)][4])):
                        same_pos_alignments[(umi,rev)] = alignment
                case _:
                    raise ValueError("Invalid option given for \"mode\"")
    else:
        #we are out of a stretch of potential duplicates
        #write out all non-duplicate reads with that value and position
        for al in same_pos_alignments.values():
            non_dupe_ctr += 1
            restore_pos(al)
            out.write("\t".join(al) + "\n")
        #reset dictionary
        same_pos_alignments.clear()
        #start tracking new position and chromosome
        same_pos_alignments[(umi,rev)] = alignment
        comparing_mol = alignment

#close files
out.close()
un_deduped.close()

print(f"Finished deduplicating file.\
        \nRemoved {bad_mapping_ctr} reads due to their being unmapped.\
        \nRemoved {bad_umi_ctr} reads due to invalid UMIs.\
        \nRemoved {duplicate_ctr} duplicate reads.\
        \nReads remaining: {non_dupe_ctr} unique reads.")
    
