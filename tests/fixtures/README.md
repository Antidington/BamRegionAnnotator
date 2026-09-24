# Regression fixture

`reads.sam` and `reference/` are synthetic, with two 1,000-base contigs and
three two-exon genes (including a reverse-strand gene). Coordinates are 1-based
in SAM/GTF and 0-based inclusive in the STAR tables. STAR chromosome starts are
0 and 1024, with a terminal offset of 2048.

`expected/*.sam` were captured from the original `src/main.rs` at
`2bbbaf75c06e4a921d0d06a54968af9f383e4607`, built using the recovered dependency
documented in [DEPENDENCIES.md](../../DEPENDENCIES.md). They were decoded with samtools 1.21.
The test suite reads these snapshots directly; it does not require samtools or
regenerate its own expected results. The four argument sets are in
`tests/baseline.rs`. This establishes wrapper regression compatibility, not
independent biological validation against a Cell Ranger release.
