#!/usr/bin/env python
#Author: André Watson, November 2023

import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description = "Deduplicate aligned single-end reads from a sorted SAM file.")
    parser.add_argument("-f", "--file", help="Absolute path to sorted SAM file", required = True)
    parser.add_argument("-o", "--outfile", help="Absolute path to output deduplicated SAM file", required = True)
    parser.add_argument("-u", "--umi", help = "File containing list of valid UMIs", required = True)
    return parser.parse_args()


def get_umi(alignment: list[str]) -> str:
    '''
    Extracts and returns an UMI from the last colon-separated element of the QNAME field of an alignment line.
    '''
    return alignment[0].split(":")[-1]

def adjust_pos(alignment: list[str]):
    '''
    Edits the value corresponding to the POS field to account for soft clipping and strandedness. Stores the old position in a valid SAM format tag.
    '''
    #Establish strandedness
    rev: bool = False
    flag: int = int(alignment[1])

    if flag & 16 == 16:
        rev = True

    #Get old position
    oldpos: int = int(alignment[3])
    newpos: int = oldpos

    #extract CIGAR string
    cigar: str = alignment[5]
    operations:list[str] = re.findall("[0-9]+[MIDNSHP=X]", cigar)

    #working with CIGAR string will be easier if we delete hard-clipped bases at the start
    if operations[0].endswith("H"):
        del operations[0]
    
    if rev:
        #sequence is being reverse complemented, adjust "rightwards"
        #ignore soft-clipping at the "start" (end of the reversed read)
        if operations[0].endswith("S"):
            del operations[0]
        for op in operations:
            if op[-1] in "MDN=XS":
                #all operations that consume reference, plus soft-clipping at the start of the reversed read
                newpos += int(op[0:-1])
        newpos -= 1 #accurate start location adjusts for off-by-one
    else:
        #sequence is not being reverse complemented, adjust "leftwards"
        #only soft-clipping at the start of the read matters for adjusting position
        if operations[0].endswith("S"):
            newpos -= int(operations[0][0:-1])

    #store old position
    alignment.append(f"OP:i:{oldpos}")
    #put new position in script
    alignment[3] = str(newpos)            


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

#define valid UMIs
umi_set: set[str] = set()
with open(umi_file,"rt") as umifile:
    for line in umifile:
        umi_set.add(line.strip())

#set up duplicate tracker and convenience variables
mols_seen: set[tuple[str,str,bool]] = set()
current_chrom: str = "-1"
bad_umi_ctr: int = 0
bad_mapping_ctr: int = 0
duplicate_ctr: int = 0
non_dupe_ctr: int = 0

#open files
out = open(outfile, "wt")
# dupes = open("./duplicates.txt", "wt")
un_deduped = open(file, "rt")

#loop over input file
while(True):
    line = un_deduped.readline().strip()
    #end of file:
    if line == "":
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
    
    #our sorted input allows us to check for duplicates by chromosome
    if alignment[2] != current_chrom:
        mols_seen.clear()
        current_chrom = alignment[2]
    #check UMI validity
    umi: str = get_umi(alignment)
    if umi not in umi_set:
        bad_umi_ctr += 1
        continue

    #get strandedness
    rev: bool = int(alignment[1]) & 16 == 16
    #adjust start position
    adjust_pos(alignment)
    if (umi, alignment[3], rev) not in mols_seen:
        #not a duplicate
        mols_seen.add((umi, alignment[3], rev))
        restore_pos(alignment)
        non_dupe_ctr += 1
        out.write("\t".join(alignment) + "\n")
    else:
        #duplicate
        # restore_pos(alignment)
        duplicate_ctr += 1
        # dupes.write("\t".join(alignment) + "\n")

out.close()
# dupes.close()
un_deduped.close()

print(f"Finished deduplicating file.\
        \nRemoved {bad_mapping_ctr} reads due to their being unmapped.\
        \nRemoved {bad_umi_ctr} reads due to invalid UMIs.\
        \nRemoved {duplicate_ctr} duplicate reads.\
        \nReads remaining: {non_dupe_ctr} unique reads.")
    
