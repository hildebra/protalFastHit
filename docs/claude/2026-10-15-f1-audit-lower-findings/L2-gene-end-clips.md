# F1 audit L2: credit for the bases clipped past a gene's end (`--gene_end_clips credit`)

2026-10-15. Branch `f1-l2-gene-end-clips`, from `f1-lower-findings` at `68316e7` (protal 0.8.5 and a doc fix); the
change is the commit that adds this report. WSL on 4 cores (`taskset -c 0-3 nice -n 5`), on a machine shared with
other sessions' builds and simulations (load 10-18).

**The finding** ([F1 audit](../2026-10-15-f1-audit-since-0.6.0a.md), L2, re-checked by a reviewer against 0.8.5):
short reads' candidates are ranked, and their MAPQ computed, by the bitscore `AlignmentInfo::Score(2, 3, 1, 2)`. A
read that runs past its gene's start or end is aligned with the gene's end fixed and the overhang free
(`AlignmentStrategy.h`, `AlignmentOrientation::Update` and the window's free ends), and `PostProcessAlignment` turns
the overhang into soft clips, which the bitscore scores 0 against 2 a match. A congener whose gene covers the whole
read therefore outranks the read's own taxon whenever its extra mismatches are fewer than 2/5 of the clipped bases:
30 bases aligned and 120 clipped score 60, a covering congener with 3 mismatches 285, MAPQ ~78 for the wrong taxon.

**The answer, in short.**
- `--gene_end_clips credit|off`, **off by default** (and for every database that records nothing): a candidate's
  bases clipped past its gene's ends, beyond those of the read's least clipped candidate, score the identity floor's
  rate (1.5 a base at `-a 0.9`) in the ranking and the MAPQ. Off, protal's outputs are byte-identical to 0.8.5's: the
  six mocks' 449,603 records, their profiles and `.truth_annotated` tables are the same as the base binary's; only the
  SAM settings line gains `gene_end_clips=off`.
- **It does what L2 asks:** with every species in the database, 6.0% fewer fragments count for a wrong species
  (5,777 → 5,431 on the 50 simulated test samples; 356 lost, 10 gained).
- **But, as designed, it dilutes the right species' MAPQs:** 2,645 fewer fragments count for the right species
  (−0.46%; 2,902 lost, mostly falling under MAPQ 4 on the right species, 257 gained), and on the training database
  (the held-out species absent) 1,816 (−0.49%) against 566 fewer wrong. The cause is the credit between candidates
  of one species, e.g. a read across two neighbouring genes of its own genome: its two partial alignments are both
  clipped at a gene end, and the more clipped one's excess credit brings it within a few points of the other
  (section 4.3; without that credit the loss is gone).
- **A per-taxon common** (every candidate of a taxon credited as that taxon's least clipped one; first a patch,
  [l2_taxon_common.patch](l2_taxon_common.patch), **since the integration the committed rule**, with its own test)
  keeps the gain and removes the loss: 340 fewer wrong
  fragments and **210 more right** ones on the finished database (39 right lost, 249 gained), 217 fewer wrong and 85
  more right on the training database; on the finished database the MAPQ of 2,843 mates falls instead of 18,179.
- **Calls** (models trained with `off`, not retrained): the same species in all six mocks; on the simulated sets the
  same 4 calls change under both rules, all at the knob (P 0.45-0.52): rows F1 0.9879 → 0.9888 (0.9885 per taxon) on
  the finished database, 0.9633 → 0.9633 (0.9629) on the training database. Within noise.
- **Recommendation:** keep `off` as the default; a database built with `credit` (the per-taxon rule) needs its models
  trained on it (the records, MAPQs and so `congener_fit_share` and the MAPQ-based features change).

## 1. The change

- `AlignmentInfo` (`src/Alignment/AlignmentUtils.h`) gains `boundary_clips`, `clip_credit` and `end_clips_credited`
  (reset in `Reset()`), and `RankScore() = Score(2, 3, 1, 2) + clip_credit`.
- `PostProcessAlignment` sets `boundary_clips`: the leading clips of an alignment that starts at gene position 0, the
  trailing clips of one that ends at the gene's length (`total_reference_length`, which `AlignAnchor` passes as
  `gene.GetLength()`; the window start it adds `del_left` to is a gene position). Clips inside a gene are none.
- `SimpleAlignmentHandler::CreditGeneEndClips` (`src/Core/AlignmentStrategy.h`), at the end of `operator()` before
  its sort, for every read and every mate: `common` = the least `boundary_clips` of the read's candidates; each
  candidate's `clip_credit = round((5 * a - 3) * (boundary_clips - common))`, the vote's floor rate (1.5 a base at
  `-a 0.9`). Equal clips earn nothing, so a read clipped alike on every candidate keeps its scores and MAPQ, which is
  what diluted the MAPQs when every clipped base was credited (the `_clip` variants of
  [the mate vote's report](../2026-10-15-mate-consensus/README.md), section 4.2).
- A mate found by `RescueMate` (the vote, or guidance) was aligned alone: `CreditFoundMate` (`src/Core/MateGuidance.h`)
  credits it against the least clipped of that mate's candidates among the fragment's.
- `RankScore()` replaces `Score(2, 3, 1, 2)` at every short-read site: the pairs' sort and MAPQ (`Bitscore` of a pair,
  `Classify.h` `SortAlignmentPairs3`, `MAPQv1`), single-end ranking and MAPQ (`ProtalSingleOutputHandler::Bitscore`),
  split mates (`BestOfMate`, `BestOfMateInTaxon`, `SettleMatesByPair`), guidance and the vote (`GuideMate`,
  `VoteMates`). With the option off `clip_credit` is 0 everywhere, so these are the old scores. Long reads keep
  `LongReadBitscore` (their segments are hard-clipped to windows), and their handler never credits.
- `AlignmentEdits` (`ZA`, and through it `congener_fit_share`) does not count a credited candidate's boundary clips as
  edits.
- **The option** `--gene_end_clips credit|off`, default `auto`: the database's `build_settings.tsv` value, else `off`.
  It is a run setting like `mate_guidance` (`BuildSettings.h`: key `gene_end_clips`, kept as text so a later
  protal's values read; `Options.cpp`: the `given` block, `FollowDatabaseRunSettings`, `RunSettingsInUse`, the help
  and the options print; `Run.cpp`: `SettingApplies` for short reads, the settings block; the SAM settings line for
  short reads, `Run/Align.cpp`; `--profile_only`'s check of a SAM's line). `--build` records the value given, else
  `off`. A database that records none was trained with `off`: the block shows `off` as the database's (no note), and
  a run given `credit` there is warned. A value this protal does not know is warned about and left at `off`.

## 2. Tests

- `tests/test_AnchoredAlignment.cpp`:
  - `GeneEndClips.PostProcessAlignmentCountsTheClipsPastTheGenesEnds`: `IIIII` + 145 `M` from gene position 0 has 5
    boundary clips, from 100 none; 145 `M` + `IIIII` to the gene's last base 5, on a longer gene none; an alignment
    whose window it does not start at (leading deletions) none.
  - `GeneEndClips.CreditedClipsRankTheReadsOwnTaxonFirst` (the finding's case, through the handler as protal calls
    it, anchored and whole-window): a random genome S; taxon 1's gene S[20..], taxon 2's S with mismatches at 40, 80,
    120; the read S[0, 150). Off: taxon 2 first, 285 against 260, the clips 20 edits. Credit: `20S130M` earns 30,
    290 against 285, taxon 1 first with MAPQ `MAPQv2(290, 285)`; its `ZA` edits 0, taxon 2's 3.
  - `GeneEndClips.EqualClipsAndReadsInsideTheirGenesEarnNothing`: both genes from S[20..] (20 clips each): no credit,
    `RankScore` = the bitscore; a read inside both genes: no clips.
  - `GeneEndClips.AFoundMateIsCreditedAgainstItsMatesCandidates`: `CreditFoundMate` off (nothing), against mate 1's
    least clips (0: 40 clips earn 60) and mate 2's (12: 42), fewer clips than the mate's candidates (nothing).
- `tests/test_RunSettings.cpp`: the key written and read back, a later protal's value kept as text and unknown, an
  empty value refused, none recorded = `off`; a run follows the database's `credit`, stays at a value given, and
  takes `off` for a database that records none or a value it does not know.
- `tests/e2e/test_protal_e2e.py` `RunSettingsTest`: the build given `--gene_end_clips credit` records it, the plain
  build records `off`; a run on it follows (`credit (auto)` in the block, `gene_end_clips=credit` in the SAM line), a
  run given `off` is warned; on the "built before" database the SAM says `off`, the block shows `off`, and a run given
  `credit` is warned.

All 626 unit tests pass (`./tests/protal_tests`), and the 191 e2e tests (`tests/e2e/test_protal_e2e.py` with
`PROTAL_TESTS_REQUIRED=1`, on a mini database built afresh by this branch's binary).

## 3. Data and method

- **Binaries:** this branch's (`--gene_end_clips off` and `credit`) and the base `68316e7`, built from `git archive`
  in WSL (Release, gcc 13).
- **Database:** the reduced GTDB r226 database `~/c4test/db/protal_db` (761 species, built by 0.8.4 at `5369c99`, pe
  gate 4 stride, models trained on `guide` alignments with the clips scored 0). It records no run settings, so every
  run here is given `--mate_guidance consensus` (protal's default for new databases). No model was retrained: the
  calls below are those of models trained with `off`.
- **Mocks:** the six paired-end mock samples in `~/Mocks/reads` (D6331 DNA, D6331 HKU and Auckland, Treichel even70,
  hCom2 in vitro, hCom1a cecum; 1M pairs each), on the finished database
  ([run_mocks.sh](../2026-10-15-mate-consensus/run_mocks.sh)).
- **Simulated test points** (the build's independent test set, `~/c4test/db/work/test/points`, 50 samples): rl150
  p10000, p100000 and p1000000 made again by [resim.sh](../2026-10-15-mate-consensus/resim.sh) (manifests and
  abundance matrices identical to the originals), rl250 p10000 and p100000 as they are; on the finished database
  (`fin`) and on the build's training database with the finished database's pe model (`tdb`: the test set's held-out
  species absent, as the build scored it) ([run_sims.sh](../2026-10-15-mate-consensus/run_sims.sh)).
- **Truth per fragment:** a read's name is `<contig>-<n>`, the contig's species from the build's work folder; a
  fragment counts for the species of its mates with MAPQ 4 or more (the profiler's cutoff)
  ([compare_records.py](../2026-10-15-mate-consensus/compare_records.py),
  [changed_fragments.py](../2026-10-15-mate-consensus/changed_fragments.py)). Calls:
  [spike_calls.py](../2026-10-15-mock-near-misses/spike_calls.py) over the `.truth_annotated` tables.
- **Scripts:** [l2_measure.sh](l2_measure.sh) runs everything, [l2_analyze.sh](l2_analyze.sh) makes the tables
  (with [l2_mapq_cross.py](l2_mapq_cross.py): the mates whose MAPQ fell or rose and those that crossed MAPQ 4); the
  per-taxon variant is this branch with [l2_taxon_common.patch](l2_taxon_common.patch) applied, run with
  `PROTAL_L2_TAXON_COMMON=1`. The tables are in [l2_results/](l2_results/) (`*_taxon.*` the variant's); the lists of
  changed fragments (`changed_*.tsv`, 0.1-1 MB) were not kept.

Commands (WSL; `L=docs/claude/2026-10-15-f1-audit-lower-findings`, `M=docs/claude/2026-10-15-mate-consensus`):

```bash
bash $M/resim.sh ~/f1l2/sims rl150_p10000 rl150_p100000 rl150_p1000000
bash $L/l2_measure.sh ~/f1l2/bin/protal ~/f1l2/basebin/protal ~/f1l2/runs
STEPS="mocks_taxon fin_taxon tdb_taxon" bash $L/l2_measure.sh - - ~/f1l2/runs ~/f1l2/expbin/protal
bash $L/l2_analyze.sh ~/f1l2/runs ~/f1l2/results
B=taxon STEPS="mocks fin tdb" bash $L/l2_analyze.sh ~/f1l2/runs ~/f1l2/results
```

## 4. Results

### 4.1 Off is 0.8.5

The six mocks with the base binary and with this branch's `--gene_end_clips off`
([identical_mocks.txt](l2_results/identical_mocks.txt)): the records of every sample (449,603, sorted) are the same,
and so are the headers apart from the settings line (which gains `gene_end_clips=off`) and `@PG`, the profiles and
the `.truth_annotated` tables. By construction too: off, no candidate gets a credit and `AlignmentEdits` counts every
clip, so `RankScore` is the old bitscore at every site.

### 4.2 Credit against off

A fragment counts for the species of its mates with MAPQ 4 or more; "right" is the species it was simulated from
([records_*.tsv](l2_results/), [mates_*.tsv](l2_results/), [changed_*.txt](l2_results/)). Credit is the committed
rule; per taxon the variant of section 4.3.

| | finished db: credit | per taxon | training db: credit | per taxon | mocks: credit | per taxon |
|---|---|---|---|---|---|---|
| fragments with records | 615,026 | | 483,780 | | 247,303 | |
| fragments whose records differ | 10,182 (1.66%) | 2,444 (0.40%) | 7,447 (1.54%) | 1,536 (0.32%) | 3,815 (1.54%) | 1,167 (0.47%) |
| of them a mate's taxon | 378 | 376 | 289 | 284 | 274 | 267 |
| crossing MAPQ 4 | 3,226 | 349 | 2,595 | 291 | 1,479 | 261 |
| mates: MAPQ down / up (same taxon) | 18,179 / 861 | 2,843 / 738 | 13,176 / 679 | 1,593 / 517 | 6,508 / 418 | 1,281 / 380 |
| mates under MAPQ 4 now / over it now (same taxon) | 5,894 / 108 | 158 / 90 | 4,704 / 117 | 168 / 69 | 2,552 / 59 | 132 / 51 |
| fragments counted for the right species (off) | 572,774 | | 370,837 | | | |
| change | **−2,645** | **+210** | **−1,816** | **+85** | | |
| fragments counted for a wrong species (off) | 5,777 | | 92,435 | | | |
| change | **−346** | **−340** | **−566** | **−217** | | |
| right lost / gained | 2,902 / 257 | 39 / 249 | 1,920 / 104 | 12 / 97 | | |
| wrong lost / gained | 356 / 10 | 350 / 10 | 690 / 124 | 322 / 105 | | |
| mates written on two taxa (off 227 / 667 / 291) | 188 | 188 | 639 | 639 | 249 | 249 |

- **The wrong fragments that go** are the finding's case. *Bifidobacterium pseudocatenulatum*
  (`rl150_p10000_s_6`, `NZ_AP012330.1-9367`): off, the mate went to *B. catenulatum*, `6S20M1X4M1X1M2X115M`, MAPQ 15;
  with credit to its own species, `36S114M`, now best but by little (MAPQ 2), so it counts for neither. *Micromonospora
  zamorensis* (`42S` and 9 mismatches on *M. sp040634005*, MAPQ 24) moves to its own species, `87S63M`, MAPQ 11, and
  counts there; *Thomasclavelia spiroformis* (`63S`, 3 mismatches and a 2-base insertion on *T. sp944324745*, MAPQ
  38) to its own, `107S43M`, MAPQ 2.
- **Most of what changes is a MAPQ, and with the committed rule mostly downwards**, on the right species: 5,894 mates
  of the finished database fall under 4 on the same taxon. Section 4.3 says why.
- On the training database credit removes more wrong fragments than the per-taxon rule (566 against 217): the held-out
  species' reads can only go to relatives, and the committed rule's lower MAPQs drop some of them under 4 as they drop
  right ones: dilution, not better placement.

**Calls** (the models were trained with `off`; [calls_*.txt](l2_results/), by the `.truth_annotated` rows):

| | off | credit | per taxon |
|---|---|---|---|
| finished db, TP / FP / FN, F1 | 1669 / 11 / 30, 0.9879 | 1672 / 12 / 26, 0.9888 | 1672 / 12 / 27, 0.9885 |
| training db | 1089 / 25 / 58, 0.9633 | 1089 / 26 / 57, 0.9633 | 1089 / 26 / 58, 0.9629 |
| mocks (six samples) | 153 / 52 / 12 | 153 / 52 / 11 | 153 / 52 / 11 |

- The finished database's 4 changed calls are the same under both rules, all at the knob: *Coprococcus sp963600665*
  (TP, P 0.484 → 0.500), *Methanococcus voltae_C* (TP, no row → 0.500), *Herpetosiphon giganteus* (TP, 0.453 →
  0.516), *Veillonella sp900549845* (FP, 0.484 → 0.500). The training database gains one FP (*Mitsuokella
  jalaludinii*, no row → 0.500).
- The mocks call the same species. hcom2 loses six rows: five absent taxa of one fragment each and one present species
  with 3 fragments (*Phocaeicola coprocola*, P 0.06, an FN either way), hence FN 12 → 11.

### 4.3 Why the committed rule dilutes, and the per-taxon common

*Megasphaera stantonii* (`rl150_p10000_s_1`, `NZ_CP029462.1-2855`; candidates with `-m 20`): mate 1 lies across the
junction of genes 2 and 117 of its own genome (neighbours in an operon), mate 2 inside gene 2. Mate 1 has two
candidates, both of the right species and both clipped at a gene end: `63S1X27M1X58M` from gene 117's first base
(164) and `46M104S` to gene 2's last (92). The least clipped is 63, so the gene-2 alignment earns 1.5 × 41 = 62 and
scores 154. The fragment's two candidates, mate 1 on gene 117 with mate 2 (164 + 295 = 459) and both mates on gene 2
(449, off 387), are then 10 apart instead of 72: MAPQ 17 → 3, and the fragment counts for nothing. The credit was meant to settle which
*taxon* fits a read that runs past a gene's end; between two places of one taxon it only adds uncertainty.

The per-taxon common ([l2_taxon_common.patch](l2_taxon_common.patch)) credits every candidate of a taxon as that
taxon's least clipped one, against the least clipped candidate of all: the taxa compete as with the committed rule,
and a taxon's own candidates keep their order and margins. In the finding's case (one candidate per taxon) both rules
are the same. On the data above it keeps 98% of the wrong fragments removed on the finished database (340 of 346)
and turns the right species' loss into a gain (−2,645 → +210; −1,816 → +85 on the training database); the mocks'
records change for 0.47% of the fragments instead of 1.54%.

## 5. Open

- **Which rule** (settled): the per-taxon common replaced the reviewer's design on the branch, with a test of its own
  (`GeneEndClips.ATaxonsCandidatesAreCreditedAsItsLeastClipped`: two candidates of one taxon at a gene junction, and a
  found mate). Its numbers are those of the `PROTAL_L2_TAXON_COMMON` runs above, the same code without the switch.
- **Models.** Every model here was trained with `off`; credit changes records, MAPQs, `ZA` and so features such as
  `congener_fit_share` and the MAPQ-based shares. A database built with `--gene_end_clips credit` trains its models on
  credited alignments, but `build_gtdb_database.py` and `collect_training_data.py` do not pass the option yet (their
  databases record `off`).
- **Scale.** Measured on the reduced r226 database (761 species) only. At GTDB scale more congeners compete, which
  makes the wrong-species case commoner, and so may the within-taxon junctions.
- **Long reads** are unchanged; their segments are clipped to windows, and their ranking (`LongReadBitscore`) was
  not part of the finding.
