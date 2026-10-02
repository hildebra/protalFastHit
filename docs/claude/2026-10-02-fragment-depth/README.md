# A fragment's bases once in the depth

- **Date**: 2026-10-02.
- **Code**: branch `audit-fixes` after `c5475ff` (with two commits of another session, `460c9a1` and `96903bd`, on
  top): `src/Profiling/Profiler.h` (`profiler::Gene`), `tests/test_StrainOutput.cpp`; docs `running.md`,
  `model-training.md`.
- **Request**: "Overlapping mates count twice in depth but once in the MSA", from
  [2026-10-02-phasing-and-foreign-genes](../2026-10-02-phasing-and-foreign-genes/README.md#which-reads-go-where).
- **Data**: the operon world's build `~/opw/b_gn2` (WSL): its 54 paired-end test samples (2x100, 2x150 and 2x250 reads
  at 1,000 to 100,000 pairs, 6 samples each), their SAMs, its training database and paired-end model.
- **Scripts**: [`fragment_depth_compare.sh`](scripts/fragment_depth_compare.sh),
  [`fragment_depth_summary.py`](scripts/fragment_depth_summary.py); runs in WSL `~/fd`.

## What changed

A taxon's depth (`Taxon::VerticalCoverage`, the abundance and the model's `depth` feature) summed, per gene, the
reference bases of every record: both mates of a pair in full, also where they overlap. The strain container counts
a fragment once (the second mate skips the part the first covered), so the MSA's coverage and the variant calls
treated two mates of one molecule as one observation, and the depth as two. Each base of the overlap was counted
twice in the depth.

`profiler::Gene::AddSam` now also keeps `m_fragment_bases`: the reference bases its fragments cover on the gene, each
once. A record of the fragment whose other mate was added to the gene before adds only what that mate did not cover
(their reference intervals' overlap is subtracted). These fragment bases, not the records' bases, now give:
- the depth from a taxon's own reads (`MappedLength`, so `VerticalCoverage` and the blended depth);
- the weights of the best reads' identity (`TopIdentity`) and the low-identity share (`LowIdentityShare`);
- a gene's depth (`Gene::VerticalCoverage`: `.profile.genes.log`'s `VCov` and `VCovExp`, the strain `.meta.tsv`'s
  `vertical_coverage`, and the features over gene depths: `depth_cv`, `stddev`, the gene variances);
- the kept conservation ratio (`conserved_fast_kept_ratio`).

`m_mapped_length` keeps every record's bases, for what counts per read: the rates per aligned kb (`lu_per_kb`,
`variant_sites_per_kb`, ...), whose numerators count both mates, `BaseIdentity`, and `.genes.log`'s
`TotalMappedLength`.

What it means: one base of depth is now a base of a sequenced molecule, as in the MSA. For a pair whose fragment F is
shorter than its two reads (2 L), the depth gets F instead of 2 L bases; per fragment at 2x150, a 220 bp fragment
gives 220/300 of its old contribution. Within a sample every taxon shares the library, so relative abundances change
only by how genes' lengths and the fragments' places on them differ between taxa. Single-end and long reads, and
pairs whose mates do not overlap, count as before.

## Effect on the test samples

Each sample profiled from its SAM (`--profile_only`) by `c5475ff` (old) and by the change (new):

| reads (fragments) | samples | TP old / new | FP old / new | FN old / new | calls changed | present taxa's depth new / old: median (5-95%) | their relative abundance new / old: median (5-95%) |
|---|---|---|---|---|---|---|---|
| 2x100 (300 +- 40) | 18 | 312 / 312 | 0 / 0 | 28 / 28 | 0 | 1.0000 (0.9980-1.0000) | 1.0002 (0.9986-1.0014) |
| 2x150 (350 +- 50) | 18 | 351 / 351 | 0 / 0 | 30 / 30 | 0 | 0.9883 (0.9682-1.0000) | 1.0032 (0.9807-1.0157) |
| 2x250 (550 +- 50) | 18 | 318 / 318 | 0 / 0 | 29 / 29 | 0 | 0.9937 (0.9761-1.0000) | 1.0017 (0.9852-1.0093) |

The collector's read setups rarely overlap (fragments a little longer than the two reads on average), so the
models trained on them see a depth 0-3% lower than before, and no call changes here. A real library of shorter
fragments (inserts of 200-300 bp at 2x150 are common) now gets a depth of its molecules' bases rather than of its
reads', so its depth feature is lower than before by up to the share of its bases in overlaps; whether that shifts
calls of low-abundance taxa there was not measured: the models are best retrained on dumps of this protal, and
`build_gtdb_database.py` does that with every build.

## Tests

Unit: `Abundance.AFragmentsOverlapCountsOnceInTheDepth` (a pair overlapping by 10 bases adds 50 of 60, a pair
apart 40, a single read 50; the records' bases still 150). Four depth tests that added separate reads under one read
id (all 0) now give each read its own: with the change their reads would have been one fragment. On the working tree
at commit time (with the other session's commits up to `07c371b`): unit 293 (292 passed, 1 skipped), end-to-end 127
passed, mini database 42 passed (2 skipped, the GTDB build tests). The comparison's table: [`results/test_samples.md`](results/test_samples.md).
