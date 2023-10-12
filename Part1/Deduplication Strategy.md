# Reference-Based PCR Duplicate Removal
## Problem Description
After alignment of sequencing reads to a reference genome, we may want to remove PCR duplicates (that is, duplicate reads arising from library amplification rather than biological processes) from our data, in order to reduce bias or incorrect interpretations in downstream analyses. This is possible with the help of Unique Molecular Identifiers (UMIs), short sequences identifying sequencing reads that originate from different molecules. 

Two sequences are PCR duplicates if all of the following are true:
+ The sequences have the same UMI
+ The sequences align to the same chromosome
+ The sequences correspond to the same strand of the reference chromosome
+ The sequences have the same 5' start position relative to the reference chromosome
  + **Note that this requirement must take strandedness and soft clipping into account (see [Initial Notes](#initial-notes)).**

## Input
We have two files as input: a list of 96 UMIs, and a sorted SAM file. `samtools sort`, by default, sorts first on the `RNAME` field (chromosome), then on the `POS` field (1-based leftmost mapping position).
## Output
We want to output one SAM file, which should be similar to the input SAM file, with the same formatting (including header lines), but having no two alignment lines that satisfy the definition of PCR duplicates. Our output SAM file should have the maximum number of alignment lines for which this condition holds.
## Proposed Algorithm
### Initial Notes
As SAM files can be prohibitively large (millions of alignment lines), we cannot simply store all the alignments in memory, contained in some data structure, and identify which alignments are PCR duplicates from there. We would prefer to develop a streaming algorithm, where we can examine one line at a time (or a sliding window of lines) and decide whether or not to retain that alignment. 

One strategy for doing so is to ensure that potential PCR duplicates are sequential in our SAM file: we can then iterate through duplicate alignment lines until we reach a non-duplicate, select and write out one of our duplicate alignments, and continue iterating (or only store and write out one duplicate, if we don't want to allow a choice of duplicates). Assuming that numbers of any single PCR duplicate are relatively low, this prevents memory issues. At first glance, `samtools sort` seems to provide our desired input. However, soft clipping and strandedness introduce several scenarios where PCR duplicates may not be apparent *or sequential*.

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
+ Open the input SAM file.
+ Open an intermediate output SAM file.
+ Read over each line of the input SAM file.
+ If the line starts with `@`, write it to the intermediate file.
+ Otherwise, we have an alignment line.
  + Split the line by tabs and store it in a list.
  + If the 0x4 bit of the `FLAG` field is set, continue on to the next line. We can't make any assumptions about unmapped reads, and we may as well check for them early.
  + If the UMI isn't in the set of known UMIs, continue on to the next line.
  + Append the value in the `POS` field to the `QNAME` field, separating it from the rest of that string with a colon. We're about to edit the POS field, so this step will let us store the original value to retrieve later.
  + If the 0x10 bit of the `FLAG` field is unset:
    + The read is not reverse complemented.
    + Check for soft clipping by examining the `CIGAR` field. For these 'forward' reads, we want to see if there are any soft-clipped bases at the start of the `CIGAR` string. Regular expression: `^([0-9]+H)?([0-9]+)S`.
    + Take the value in the `POS` field. Subtract the amount of soft-clipped bases at the start of the `CIGAR` field, if any, and store that value back in the `POS` field.
    + Write the line to the intermediate file.
    + Continue to the next line.
  + Otherwise, if the 0x10 bit of the `FLAG` field is set:
    + The read is being reverse complemented.
    + Check for soft clipping by examining the `CIGAR` field. We actually still only need the amount of soft-clipped bases at the start of the `CIGAR` string.
    + Take the value in the `POS` field. Subtract the amount of soft-clipped bases at the start of the `CIGAR` field. Add the length of the string in the `SEQ` field. Subtract 1. Store that value back in the `POS` field.
    + Write the line to the intermediate file.
    + Continue to the next line.
+ When there are no lines left to read in the input SAM file, close both files.
### Phase 2: Resorting
+ Run `samtools sort` on our intermediate file. Specify an output SAM file.
+ We haven't violated any formatting rules (just made the content temporarily inaccurate), so this should be possible. Potential PCR duplicates should now be next to each other. 
### Phase 3: Filtering
#### 3a: UMI Storage
+ Create a set to store UMIs.
+ Open the input list of UMIs.
+ Read each line of the file.
+ Add the contents of each line to the set.
+ Close the file with the list of UMIs.
#### 3b: Writing Output
+ Create an empty list to store an alignment line.
+ Open the output from Phase 2.
+ Open a final output SAM file.
+ Read each line of the Phase 2 file.
+ If the line starts with `@`, write it to the final output SAM file.
+ Otherwise, we have an alignment line.
+ Split the line into a list by tabs.
+ If this is the first alignment line we see:
  + Store it in the list we made earlier.
  + Continue to the next line.
+ If this is not the first alignment line:
  + If all of the following are true:
    + This line and the stored line have the same UMI
    + This line and the stored line have the same value in the `POS` field
    + This line and the stored line have the same value in the `RNAME` field
    + The 0x10 bits in this line's `FLAG` and the stored line's `FLAG` are either both set or both unset
  + Do not store this line. Continue to the next line.
  + If any of those conditions were not met:
    + Replace the value in the stored line's `POS` field with the value we appended to the `QNAME`.
    + Remove that value (and the colon we added) from the `QNAME` field of the stored line.
    + Write the stored alignment line to our final output file.
    + Replace the stored line with the current line.
    + Continue to the next line.
+ If there are no more lines in the file:
  + Replace the value in the stored line's `POS` field with the value we appended to the `QNAME`.
  + Remove that value (and the colon we added) from the `QNAME` field of the stored line.
  + Write the stored alignment line to our final output file.
+ Close all open files.
#### 3c: Cleanup
+ Delete all intermediate files.

## Necessary Helper Functions
```python
def validate_umi(alignment: list) -> bool:
  '''
  Checks to see if an (unedited) alignment line has an UMI that appears in our list of valid UMIs. Returns True if the UMI is valid, False otherwise.
  '''
  Input Example:
  ["NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC", "0", "2", "76814284", "36", "71M", "*", "0", "0", "TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA", "6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/", "MD:Z:71", "NH:i:1", "HI:i:1", "NM:i:0", "SM:i:36", "XQ:i:40", "X2:i:0", "XO:Z:UU"]
  Output Example: True

```


```python
def adjust_pos(alignment: list):
  '''
  Edits the value corresponding to the POS field to account for soft clipping. Stores the old position in the first field (QNAME).
  '''
  Input Example:  ["NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC", "0", "2", "76814284", "36", "2S69M", "*", "0", "0", "TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA", "6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/", "MD:Z:71", "NH:i:1", "HI:i:1", "NM:i:0", "SM:i:36", "XQ:i:40", "X2:i:0", "XO:Z:UU"]
  End State Example: ["NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC:76814284", "0", "2", "76814282", "36", "2S69M", "*", "0", "0", "TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA", "6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/", "MD:Z:71", "NH:i:1", "HI:i:1", "NM:i:0", "SM:i:36", "XQ:i:40", "X2:i:0", "XO:Z:UU"]
```

```python
def restore_pos(alignment: list):
  '''
  The reciprocal operation of adjust_pos(). Takes an edited alignment line, strips the original POS value off of the first field, and uses it to replace the edited POS value.
  '''
  Input Example:  ["NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC:76814284", "0", "2", "76814282", "36", "2S69M", "*", "0", "0", "TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA", "6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/", "MD:Z:71", "NH:i:1", "HI:i:1", "NM:i:0", "SM:i:36", "XQ:i:40", "X2:i:0", "XO:Z:UU"]
  End State Example: ["NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC", "0", "2", "76814284", "36", "2S69M", "*", "0", "0", "TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA", "6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/", "MD:Z:71", "NH:i:1", "HI:i:1", "NM:i:0", "SM:i:36", "XQ:i:40", "X2:i:0", "XO:Z:UU"]

```

```python
def detect_duplicate(a1: list, a2:list) -> bool:
  '''
  Returns True if the two (adjusted) alignments are duplicates, False if they aren't.
  '''
  Input Example:
  Output Example: True

```
