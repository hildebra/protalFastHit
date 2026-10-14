# Why a grid-picked sparse index calls fewer strains at C 4 when its reads align the same

The user asked why paired-end calls drop on a sparse index (an index holding 1/C of the cores, not only a gate looking
up 1/C of them), when the grid's C 4 runs report about the same reads; and whether at C 4 on the sparse index the number
of aligned reads differs much. The case is follow-up 5 of
[2026-10-13-gene-grid-index](../2026-10-13-gene-grid-index.md): on the strain world (532 held-out real strains), the
grid at C 4 on its sparse index (`g4`) wrote the same records as the ungated run (`g0`), but the shared presence model
called 0.867 of the test strains against 0.900.

**Answer.** The aligned reads do not differ: 176,940 primary records against 176,958 (−0.01%), the same share aligned
and on the right species. What differs is one per-read count, the read's seeds on unique k-mers (SAM tags `ZU`, `ZT`).
A read seeded from its grid cores alone has a quarter of its seeds, so it counts a quarter of its unique k-mers; the
features made from them (`lu_per_kb`, `lsu_per_kb`, `lu_gene_rate3`, `lsu_gene_rate3`) fall to 0.26-0.75 of the ungated
run's, and a model trained on fully seeded runs reads that as a relative's reads. Putting the ungated run's `ZU`/`ZT`
into the sparse run's SAMs brings back every lost call (0.907); multiplying the sparse run's own counts by 4 does too
(0.900, the same false positives). Every other difference (failed candidates tripled by the single-seed anchors, MAPQ,
crowding, untried candidates, alternatives) changes no call of this model.

## Data and commands

- Commit `47e3cc2` (built from `git archive` in WSL `~/grid-verify`, the binary kept as `~/strainfn/cf/protal`). The
  profiling code is that of the runs' binary (`19ccc0e` + the whole-index change) for this model's features: the two
  runs re-profiled unchanged give exactly their calls (40 of 40 samples).
- The strain world of [2026-10-11-strain-false-negatives](../2026-10-11-strain-false-negatives/README.md) in WSL
  `~/strainfn`: database `db_grid`, 20 deep and 20 shallow paired-end samples, the small real world's paired-end GBM
  (`model_s1/model_pe.xml`, trained on ungated, fully seeded runs). The SAMs are those of follow-up 5's runs `g0`
  (`--sparse_first 0`) and `g4` (`--sparse_first 4 --sparse_first_pick grid`, every sample grid-picked: the sparse
  index, single-seed anchors).
- `samdiff_all.sh` (`samdiff.pl` per sample): the two runs' SAMs record by record (read name and flag).
- `cf_runs.sh name:base:donor:spec ...`: `swap_tags.pl` writes the base run's SAMs with the fields named in `spec` taken
  from the donor run's same record (or `ZU`/`ZT` multiplied, `x4`); `protal --profile_only` with `--profile_truth` (the
  samples' species from `reads/truth.tsv`) profiles them with the same model and flags as the runs (`--no_strains
  --no_qcmsa`), writing the training dumps. Run: `cf_runs.sh g0:g0:g0:none g4:g4:g4:none`, then `zu:g4:g0:ZU,ZT
  zf:g4:g0:ZF zuzf:g4:g0:ZU,ZT,ZF zn:g4:g0:ZN mapq:g4:g0:MAPQ zc:g4:g0:ZC za:g4:g0:ZA
  all:g4:g0:MAPQ,ZU,ZT,ZA,ZF,ZC,ZN,ZR x4:g4:g4:x4 x4zf:g4:g0:x4,ZF`.
- Scored with that study's `evaluate.py ~/strainfn ~/strainfn/eval_cf cf_g0 cf_g4 cf_zu ...` (calls of the shallow
  samples, false positives of shallow and deep); features compared with `compare_dumps.pl run_cf_g0/profiles
  run_cf_g4/profiles model_s1/model_pe.xml`.

## The reads: the same

All 40 samples, `g0` against `g4`:

| | `g0` | `g4` |
|---|---|---|
| primary records with an alignment | 176,958 | 176,940 (−0.01%) |
| pairs given up by the gate (follow-up 5) | 0 | 131 of 88,800 |
| candidate alignments (follow-up 5) | 452,927 | 557,934 (+23%) |
| records on another gene (of 176,829 in both) | | 1,218 (0.7%) |
| records with another MAPQ | | 31,856 (18%; 18,390 lower, 13,466 higher) |
| of them across MAPQ 4 (the profiler's filter) | | 487 (0.28%) |
| `ZU` summed (seeds on unique k-mers) | 2,126,401 | 553,348 (0.260) |
| `ZT` summed (unique at distance two) | 1,749,927 | 520,334 (0.297) |
| `ZN` summed (taxa the seeds cannot tell apart) | 286,840 | 272,233 (0.949) |
| records with another `ZF` (failed candidates) | | 48,232 (27%) |

The read level of follow-up 5 says the same: aligned 0.996 / 0.999 / 0.984 (easy / medium / hard) in both, right
species 0.907 against 0.897 for the hard strains (`g4` ahead). The extra candidates are failed alignments of the
single-seed anchors, not more reads.

## The features: the unique k-mer counts at a quarter

True species with a profile row in both runs (1,777); the model's features that changed in more than 1% of them,
the largest shown:

| feature | changed | mean `g0` | mean `g4` | median `g4`/`g0` | the 31 lost: `g0` | `g4` |
|---|---|---|---|---|---|---|
| `lu_per_kb` | 99% | 81.1 | 21.2 | 0.260 | 40.7 | 9.9 |
| `lsu_per_kb` | 98% | 66.0 | 19.9 | 0.314 | 28.2 | 9.0 |
| `failed_candidate_rate` | 92% | 0.045 | 0.194 | 3.20 | 0.218 | 0.541 |
| `lu_gene_rate3` (genes with > 5 unique seeds) | 83% | 0.259 | 0.198 | 0.750 | 0.057 | 0.015 |
| `failed_gene_share` | 83% | 0.039 | 0.162 | 2.33 | 0.211 | 0.511 |
| `lsu_gene_rate3` | 83% | 0.249 | 0.192 | 0.750 | 0.043 | 0.012 |
| `mean_mapq` | 73% | 56.0 | 54.8 | 0.991 | 40.9 | 41.9 |
| `seed_crowding` | 72% | 0.369 | 0.316 | 0.869 | 0.649 | 0.562 |

The others move by a few percent at most (`lu_gene_rate2`, `lsu_gene_rate2`, `low_mapq_share` and the congener-fit
shares by 1-4%, `untried_candidate_rate` 0.022 to 0.019); depth, fragments, identity, coverage and the allele and
ancestry features by less than 1%. The model's probability of the true species falls from 0.969 to 0.949 on average,
and for the 31 it loses (shallow and deep samples; 6 gained) from 0.768 to 0.306.

## Which change costs the calls: one tag at a time

`g4`'s SAMs with one field taken from `g0`'s same records, re-profiled (shallow samples; 600 strain instances, by their
pairs; FP: called species absent from the sample):

| SAMs | 1-2 pairs | 3-6 | 8-16 | 24-32 | all | hard | FP shallow / deep |
|---|---|---|---|---|---|---|---|
| `g0` | 0.783 | 0.894 | 0.939 | 0.967 | 0.900 | 0.842 | 2 / 12 |
| `g4` | 0.675 | 0.878 | 0.928 | 0.950 | 0.867 | 0.743 | 1 / 7 |
| `g4`, `g0`'s `ZU`, `ZT` | 0.783 | 0.911 | 0.944 | 0.967 | **0.907** | 0.851 | 2 / 12 |
| `g4`, `g0`'s `ZF` | 0.683 | 0.861 | 0.928 | 0.950 | 0.863 | 0.723 | 1 / 7 |
| `g4`, `g0`'s `ZU`, `ZT`, `ZF` | 0.800 | 0.906 | 0.944 | 0.967 | 0.908 | 0.861 | 2 / 12 |
| `g4`, `g0`'s MAPQ | 0.675 | 0.861 | 0.928 | 0.958 | 0.863 | 0.723 | 1 / 8 |
| `g4`, `g0`'s `ZN` / `ZC` / `ZA` | 0.675 | 0.878 | 0.928 | 0.950 | 0.867 | 0.743 | 1 / 7 |
| `g4`, all of `g0`'s tags and MAPQ | 0.783 | 0.889 | 0.939 | 0.967 | 0.898 | 0.832 | 2 / 12 |
| `g4`, its own `ZU`, `ZT` × 4 | 0.758 | 0.911 | 0.939 | 0.967 | **0.900** | 0.842 | 2 / 12 |
| `g4`, its own `ZU`, `ZT` × 4, `g0`'s `ZF` | 0.792 | 0.900 | 0.939 | 0.967 | 0.903 | 0.861 | 2 / 12 |

- **The unique k-mer counts are the whole loss.** With `g0`'s counts the sparse run calls 0.907 (all classes) and
  0.851 (hard), at `g0`'s false positives; the misses the model turned down ("2+ pairs counted, model said no") go from
  27 to 13 (`g0` 17), those with one pair counted from 30 to 20 (`g0` 20).
- **Scaling them by C is enough here**: `ZU`, `ZT` × 4 calls 0.900 / 0.842 with 2 / 12 false positives, `g0`'s figures
  exactly (other strains: 0.758 against 0.783 for 1-2 pairs, ahead at 3-6). At lower knobs the scaled run is within
  0.01 of `g0` (knob 0.2: 0.937 against 0.940, 6 against 8 false positives).
- **The tripled failed-candidate rate costs nothing with this model**: `g0`'s `ZF` alone changes 0.867 to 0.863. A
  model that weighs `failed_candidate_rate` more (the r226 models use it against relatives' calls) may react to it.
- MAPQ, crowding, untried candidates and alternatives change at most a call or two.
- The fewer false positives of `g4` (1 / 7 against 2 / 12) come from the same counts: with `g0`'s or scaled ones they
  are `g0`'s again. So the sparse index did not make the evidence cleaner; the model is stricter with every species.

## Why the counts scale and the reads do not

A read's record depends on its anchors: once a gene is anchored the whole mate is aligned against it, and since
follow-up 3 one grid seed anchors a short read's gene. So the records are the same with a quarter of the seeds. `ZU`
and `ZT` are counted from the seeds of the read's anchor (`ChainAnchorFinder::CountUniques`), which on a sparse index
are the read's grid cores alone: about 1/C of its unique k-mers (0.26 here; `ZT` 0.30, probably because the grid prefers
cores of low key ubiquity). The reference's side of these features (`lu_rate_ref`, the
genes with unique k-mers that `lu_gene_rate*` divides by) is counted from the whole reference and does not change, so
per kb and per gene the read side comes out at a quarter. The same was seen with the hash pick's emulated 1/8 index
on HiFi ([2026-10-13-sparse-seeding](../2026-10-13-sparse-seeding/README.md): `lu`/`lsu` features / 8, 8 true positives
lost with the same reads).

Where reads are lost, the calls drop for that reason too: on the mini database the hash pick's emulated 1/4 index lost
6.5% of the paired-end records and the grid with two-seed anchors 3.2% (follow-ups 2 and 3 of the grid report); with
single-seed anchors the grid loses 0.05%, and here none.

## What it means

- A model trained on fully seeded runs cannot score a sparse-index run as it is: its unique k-mer features are at
  1/C. Either the counts are made independent of the seeding, or the models are trained on runs seeded as the sparse
  index (`build_gtdb_database.py --sparse-first-pick grid`; its collections seed so, `--sparse_index emulate`).
- Making them independent is simple: on a sparse index (and with `--sparse_index emulate` or `--sparse_seeding`), write
  `ZU` and `ZT` multiplied by the index's ratio of unique entries in the whole index to those of the seed level
  (about C; the sparse load decodes every entry, so it can count both). That restored this run exactly; the counts of a
  read become a 1/C sample scaled up, noisier for a species of one or two reads. The r226 `--test-grid` rows, scored by
  models trained otherwise, would then measure the grid rather than the model's mismatch.
- The failed candidates of the single-seed anchors (`failed_candidate_rate` × 3) did nothing here but should be
  checked with an r226 model; leaving single-seed anchors out of `ZF` would keep that feature as it is in other runs.
- Not built; this report only measures.

## Follow-up: which features depend on read length or on the seeding (2026-10-14)

The user proposed to drop k-mer counts from the models, as they tie a model to a read length, and to use instead the
fraction of the tested cores (potential anchors) that hit the species, with `ZU` and `ZT` as fractions and the
number of tested cores as features; then asked which other features depend strongly on read length or on the number
of anchors, and suggested the aligned read length as a feature.

**Data and commands.** As above (`47e3cc2`, `~/strainfn/cf/protal`, `db_grid`, `model_s1`). Read length:
`len_runs.sh 150 100 75` cuts every mate of the strain world to its first L bases (the same fragments, positions and
read names; the same pairs, so the depth falls with L; the reads' simulated errors grow along the read, so the cut
reads also have fewer errors) and runs them ungated with training dumps. `shift_all.sh` (`feature_shift.pl`) compares
every column of the dumps, 150 against 100 and 75 bp and `g0` against `g4`: shift = the change of the mean over the
present (or absent) species with a row in both runs, divided by the gap between present and absent species in the first
run (a shift of 1 moves a feature as far as presence does; features with almost no gap look large on this scale); and
the median ratio over the present species. `cf_runs.sh` again for tags of the 150 bp run's same read in the 75 bp run's
SAMs (`same`: MAPQ only from a record on the same gene, 95.5% of them), `evaluate.py` for the calls, `extra_fp.pl` for
the false positives the 75 bp run adds.

### The denominator of a unique-seed fraction

Summed over the 40 samples (the runs' seeding lines; `ZU` from the SAMs):

| | `g0` (ungated) | `g4` (sparse index, C 4) | `g4` / `g0` |
|---|---|---|---|
| cores looked up | 3,977,957 | 3,972,220 | 1.00 |
| of them in the index | 3,693,757 | 1,428,396 | 0.39 |
| seeds (all taxa) | 7,161,266 | 2,363,204 | 0.33 |
| `ZU` | 2,126,401 | 553,348 | 0.26 |
| `ZU` per looked-up core / per core in the index / per seed | 0.535 / 0.576 / 0.297 | 0.139 / 0.387 / 0.234 | 0.26 / 0.67 / 0.79 |

On the sparse index every core is still looked up, so the cores tested are no denominator: per tested core the
unique count is at 0.26 as before. Per seed it is at 0.79: the grid keeps conserved cores (they find more seeds,
0.33 rather than 0.25) that are less often unique. A fraction of the seeds needs the reference's unique share among
the index's loaded entries as its expectation, to come out the same on the sparse index.

### Calls by read length (the 150 bp model)

| reads | strains called | hard | FP shallow / deep | hard strains' pairs: right species / counted on another |
|---|---|---|---|---|
| 150 bp | 0.900 | 0.842 | 2 / 12 | 0.897 / 0.038 |
| 100 bp | 0.903 | 0.891 | 9 / 21 | 0.868 / 0.060 |
| 75 bp | 0.915 | 0.871 | 13 / 39 | 0.825 / 0.088 |

Shorter reads put more of a strain's pairs on another species, and the model, trained on 150 bp reads, calls more of
those species: the opposite error to the sparse index's.

### Features that move

The model's features (and those of the default set) whose shift exceeds 0.1 of their present/absent gap, or whose
present species' median ratio moves by more than 10%, grouped by cause (ratios: present species' median; arrows:
absent species' means):

| cause | features | 150 -> 100 bp | 150 -> 75 bp | `g0` -> `g4` (seeding) |
|---|---|---|---|---|
| seeds on unique k-mers | `lu_per_kb`, `lsu_per_kb` | about 1 | 0.79, 0.79 | 0.26, 0.31 |
| | `lu_gene_rate3`, `lsu_gene_rate3` | 0.92 (lsu) | 0.87, 0.76 | 0.75, 0.75 |
| candidates aligned | `failed_candidate_rate`, `failed_gene_share` | 0.34 -> 0.23 | 0.34 -> 0.15 | 0.34 -> 0.65 (present × 3.2) |
| | `untried_candidate_rate` | little | 0.21 -> 0.15 | 0.19 -> 0.14 |
| MAPQ | `mean_mapq` | 36 -> 43 | 35 -> 47 | unchanged |
| | `uniqueness` (share of records above MAPQ 20) | 0.46 -> 0.65 | 0.46 -> 0.77 (present at 150 bp: 0.57) | unchanged |
| | `low_mapq_share` | 0.57 -> 0.49 | 0.58 -> 0.43 | unchanged |
| the best read | `top_identity`, `allele_top_identity` | 0.976 -> 0.982 | 0.976 -> 0.986 | unchanged |
| ties with congeners | `congener_fit_share` | 1.52 | 2.00 | unchanged |
| | `seed_crowding` | 1.13 | 1.18 | 0.87 |
| | `split_fragment_share`, `mate_lost_share` | 0.029 -> 0.070, 0.060 -> 0.074 | 0.027 -> 0.102, 0.052 -> 0.093 | unchanged |
| per record | `ancestry_sites_per_record`, `ancestry_indel_sites_per_record` | 0.66, 0.67 | 0.49, 0.40 | unchanged |
| depth per position | `variant_sites_per_kb` | 0.52 | 0.28 | unchanged |
| | `multiallelic_sites_per_kb`, `fixed_difference_rate`, `polymorphic_site_rate` | about half (polymorphic: 0) | 0 (none left) | unchanged |
| | `breadth_ratio`, `gene_divergence_dispersion`, `RAF1` | 0.97, 0.82, 0.62 | 0.95, 0.74, 0.37 | unchanged |
| depth | `depth` | 0.67 | 0.50 | unchanged |

Little or nothing moves in the mean identity, coverage-shape, conservation, ancestry-agreement, allele, column-weight
and relative (`em_*`, `relative_*`) features. The counts outside the default set (`lu`, `lsu`, `lu_per_read`,
`unique_hits`, `A*`, `AF*`) move most of all.

### Which of them make the 75 bp false positives

The 75 bp run's SAMs with fields of the 150 bp run's same read, re-profiled:

| SAMs | called | FP shallow / deep |
|---|---|---|
| 150 bp | 0.900 | 2 / 12 |
| 75 bp | 0.915 | 13 / 39 |
| 75 bp, MAPQ of 150 bp (same gene) | 0.912 | 7 / 34 |
| 75 bp, `ZA` / `ZF` / `ZF`+`ZC` of 150 bp | 0.917 | 13 / 39-40 |
| 75 bp, MAPQ, `ZA`, `ZF`, `ZC`, `ZN`, `ZR` of 150 bp | 0.917 | 7 / 34 |

- **MAPQ makes half of the shallow extras.** A short read seeds fewer anchors on congeners, so fewer alternatives are
  aligned and its MAPQ is higher although it tells the species apart less well.
- **The rest are reads that now look perfect.** All 31 species the 75 bp run (with 150 bp MAPQ) adds have a best read
  of identity 1. At 150 bp, 16 of them had 0.973-0.993, 7 had 1 already (and were turned down on other features), 8 had
  no record. A relative's read at 99.3% identity is free of differences with probability 0.993^150 = 0.35 at 150 bp and
  0.59 at 75 bp. Their fragments and mean identity hardly differ; 23 of the 31 are *Priestia megaterium* (and `_C`) or
  *Dorea_A longicatena* (and `_C`), congeners of the strain study's two species that are near-identical to a congener
  (*P. megaterium_D*, *D. amylophila*). `top_identity` is a maximum over reads, so
  it grows with fewer bases per read (and with more reads).
- Taking the 150 bp `ZU`/`ZT` as well raised the false positives to 20 / 76 (counts of 150 bp seeds per kb of a 75 bp
  alignment, twice too high): not a valid swap, left out.

### What it means for the plan

- **Fractions instead of k-mer counts: yes**, with the anchor's seeds as the denominator (a new tag: the seeds of the
  record's anchor) rather than the cores tested, and the reference's unique share among the loaded entries of the
  record's genes as its expectation. `lu_per_kb` also loses 21% at 75 bp because a seed needs its flanks (fewer seeds
  per base in a short read); a share of the seeds does not.
- **The number of tested cores as a feature: better not as a count.** It moves with read length, C and error rate, which
  is what the fractions remove; it says how much a fraction is worth, which the taxon's fragments and the read length
  (below) say as well.
- **Other features by construction:** sites per aligned kb instead of per record (ancestry); the best read's identity
  as a lower bound given its aligned length (a perfect 75 bp read ~0.96, a perfect 150 bp read ~0.98, by the rule of
  three), so that a short perfect read does not read as 1; failed candidates without the single-seed anchors' (the
  seeding part of `failed_candidate_rate`).
- **The aligned read length as a feature: yes, for what cannot be made invariant.** MAPQ, ties with congeners, the
  candidates a read finds and split pairs carry less information in a short read, and only a model that knows the
  length can weigh them. As a sample feature (the median aligned bases of the sample's records, a library property as
  `sample_log_fragments` is the depth), and the taxon's own against it (partial alignments at gene ends, clipped
  relatives). It only helps if the training samples span read lengths (collections drawing each sample's read length);
  trained on 150 bp alone the model never sees it vary.
- Not built. The depth rows come from keeping the pairs while cutting the reads; the training's varied depths already
  cover depth.
