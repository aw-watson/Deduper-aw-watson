#!/usr/bin/env python
#Author: André Watson, November 2023

import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(description = "Adjust start positions in a SAM file to account for strandedness and soft-clipping")
    parser.add_argument("-f", "--file", help="Absolute path to input SAM file", required = True)
    parser.add_argument("-o", "--outfile", help="Absolute path to output adjusted SAM file", required = True)
    return parser.parse_args()


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


#main script

#read in arguments
args = get_args()
file : str = args.file
outfile: str = args.outfile

#open files
out = open(outfile, "wt")
un_adjusted = open(file, "rt")

while(True):
    line = un_adjusted.readline().strip()
    if line == "":
        break
    if line.startswith("@"):
        out.write(line + "\n")
        continue
    alignment = line.split("\t")
    adjust_pos(alignment)
    out.write("\t".join(alignment) + "\n")

out.close()
un_adjusted.close()


