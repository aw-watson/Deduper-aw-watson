# Reference-Based PCR Duplicate Removal
## Problem Description
After alignment of sequencing reads to a reference genome, we may want to remove PCR duplicates (that is, duplicate reads arising from library amplification rather than biological processes) from our data, in order to reduce bias or incorrect interpretations in downstream analyses. This is possible with the help of Unique Molecular Identifiers (UMIs), short sequences identifying sequencing reads that originate from different molecules. 

Two sequences are PCR duplicates if all of the following are true:
+ The sequences have the same UMI
+ The sequences align to the same chromosome
+ The sequences correspond to the same strand of the reference chromosome
+ The sequences have the same 5' start position relative to the reference chromosome
  + **Note that this requirement must take strandedness and soft clipping into account (see [Initial Notes]()).**

## Input
We have two files as input: a list of 96 UMIs, and a sorted SAM file. `samtools sort`, by default, sorts first on the `RNAME` field (chromosome), then on the `POS` field (1-based leftmost mapping position).
## Output
We want to output one SAM file, which should be similar to the input SAM file, with the same formatting (including header lines), but having no two alignment lines that satisfy the definition of PCR duplicates. Our output SAM file should have the maximum number of alignment lines for which this condition holds.
## Proposed Algorithm
### Initial Notes
As SAM files can be prohibitively large (millions of alignment lines), we cannot simply store all the alignments in memory, contained in some data structure, and identify which alignments are PCR duplicates from there. We would prefer to develop a streaming algorithm, where we can examine one line at a time (or a sliding window of lines) and decide whether or not to retain that alignment. 

One strategy for doing so is to ensure that potential PCR duplicates are sequential in our SAM file: we can then iterate through duplicate alignment lines until we reach a non-duplicate, select and write out one of our duplicate alignments, and continue iterating. Assuming that numbers of any single PCR replicate are relatively low, this prevents memory issues. At first glance, `samtools sort` seems to provide our desired input. However, soft clipping and strandedness introduce several scenarios where PCR duplicates may not be apparent *or sequential*.

Soft clipping is a feature of the SAM format specification. The `POS` field gives the 1-based leftmost mapping position for the alignment, but this is not necessarily the starting position of the sequencing read, only the starting position of the aligned bases relative to the 5' end of the reference chromosome. At the beginning and end of the string in the `CIGAR` field, we can recover encoded counts for how many bases did not align from the start or the end of the read. Soft-clipped bases are not accounted for by `POS`.

Consider alignment lines (ALs) with the following (abbreviated) fields:
+ ...
+ AL 1: `RNAME:1  POS:10  CIGAR:1S9M`
+ AL 2: `RNAME:1  POS:11  CIGAR:2S7M1S`
+ AL 3: `RNAME:1  POS:50  CIGAR:1S9M`
+ AL 4: `RNAME:1  POS:52  CIGAR:10M`
+ AL 5: `RNAME:1  POS:54  CIGAR:5S5M`
+ ...

Lines 1 and 2 are potentially duplicates, but that is not immediately apparent. If we account for soft-clipped bases, we see that had the reads fully aligned, they would have both started at base position 9 on the reference chromosome. 
Lines 3 and 5 are also potentially duplicates, with the first base of the reads corresponding to base position 49 on the reference chromosome. However, sorting on the `POS` field will render them non-sequential.

Strandedness complicates our strategy further. If a read aligns to the complementary strand of the reference chromosome, the 0x10 bit (2<sup>4</sup>) in the `FLAG` field will be set. The read is reverse complemented in our SAM file, so its real start position is the end of the sequence (closer to the 3' end of the reference chromosome), which the `POS` does not accurately correspond to. Combined with soft clipping, this can again cause issues with detecting duplicates.

Consider alignment lines with the following (abbreviated) fields:
+ ...
+ AL 1: `FLAG:000000010000  RNAME:1  POS:30  CIGAR:8M3S`
+ AL 2: `FLAG:000000010000  RNAME:1  POS:35  CIGAR:1S2M3S`
+ ...

Lines 1 and 2 are potential duplicates! The sequencing reads would have started in the same place (base position 40 on the reference chromosome) when both strandedness and soft clipping are taken into account. This wouldn't be as much of an issue if sequences were always the same length, but that's not guaranteed.

The first two phases of this proposed process account for these cases to create an intermediate file with a more ideal ordering of alignment lines; the third phase processes this intermediate file into our final output.

### Phase 1: Adjusting Input
### Phase 2: Resorting
### Phase 3: Filtering
#### 3a: UMI Storage
#### 3b: Writing Output
#### 3c: Cleanup

## Necessary Helper Functions
