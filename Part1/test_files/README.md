## Guide to test files:
UND\[#\]: UMI non-duplicates: different UMIs mean these reads are not duplicates.

SND\[#\]: Strand non-duplicates: different strands mean these reads are not duplicates.

LND\[#\]: Location non-duplicates: different starting positions mean these reads are not duplicates.

RND\[#\]: RNAME non-duplicates: these reads are on different chromosomes and are not duplicates.

SD\[#\]: simple duplicates: these reads are easily identifiable duplicates.

SCD\[#\]: soft-clipped duplicates: when adjusted for soft clipping, these reads are duplicates.

NSD\[#\]: non-sequential duplicates: when sorting by `POS`, reads NSD1 and NSD3 are not placed next to each other, but adjusting for soft clipping reveals them to be duplicates.

RD\[#\]: reversed duplicates: both of these reads are reverse complemented, and would have started at the same position on the complementary strand of the reference sequence. They are duplicates.

RSCD\[#\]: reversed, soft-clipped duplicates: As with the previous set of duplicates, but with soft-clipped bases to take into account.

RVND\[#\]: reversed, naive duplicate, non-duplicates: These reads appear to be duplicates when looking at `POS`, but are not actually duplicates, since they start at different positions when accounting for them both being reverse complemented.

RSCVND\[#\]: reversed, soft-clipped, naive duplicate, non-duplicates: As with the previous set of non-duplicates, but with soft-clipped bases to take into account.
