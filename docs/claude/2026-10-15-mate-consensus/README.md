# Mate guidance as a vote of both mates (`--mate_guidance consensus`)

2026-10-15. Branch `worktree-agent-a7bdc6cbd67cec5a5`, from `audit-fixes` at `dd01eb1`; the change is the commit that
adds this report. WSL on 4 cores (`taskset -c 0-3 nice -n 5`), on a machine shared with other sessions' builds and
benchmarks (load 5-22).

**The question** (the user's): "revisit the mate guidance: the best situation would be a consensus vote from both
mates, that is weigh the evidence." It follows finding L1 of the [F1 audit](2026-10-15-f1-audit-since-0.6.0a.md):
guidance completed a pair on the sure mate's taxon without weighing the other mate's own best alignment.

**The answer, in short.**
- `--mate_guidance consensus|guide|off`, default `consensus`. `guide` is protal 0.7.1-0.8.4's guidance, unchanged
  (its records are those of the build's own 0.8.4 runs, 447,059 fragments compared); `off` is `--no_mate_guidance`,
  which stays as its alias.
- The vote changes 0.9-2.0% of the fragments' records. With the models trained on `guide` alignments the calls stay
  the same within 1-2 species of ~2,080 present on each simulated test set (50 samples), and do not change at all
  on the six paired-end mocks. Fewer fragments count for the wrong species: −0.3% where every species is in the database, −0.4% where
  the test set's held-out species are absent. Slightly fewer count for the right species (−0.02% to −0.05%), mostly
  fragments that drop from MAPQ 4 to 3.
- On the build's own independent test set, scored as the build scores it, F1 goes from 0.9597 (TP 678, FP 18, FN 39:
  exactly the build's figures) to 0.9618 (680, 17, 37). That is within noise.
- Cost: +2.65% of the paired-end alignment loop's instructions against `guide` (+3.7% against `off`) on a real sample (200,000 pairs of D6331 HKU).
- **Recommendation: `consensus` as the default** (set in this commit). It is neutral to slightly positive on every
  measure here, and it removes the confident wrong pairs that L1 describes. The databases' models need no retraining
  for it: their features move little and only taxa at the knob's edge change call.
- The setting must be recorded in the database's `build_settings.tsv` (another session adds a `mate_guidance` key,
  a string with exactly these values `consensus|guide|off`), so that a run can tell whether its alignments match the
  ones the models were trained on. A database built before has no key; its models were trained with `guide`.

## 1. What `guide` did (L1, checked in the code)

At `dd01eb1`, after a pair's candidates were joined and sorted by the sum of their mates' scores (Bitscore), a best
candidate with one mate only went to `GuideMate` (`src/Core/MateGuidance.h`):
- the mate with an alignment guided when its MAPQ against **its own** other alignments was 20 or more and the other
  mate had **no candidate of its taxon** (`other_has_taxon`, lines 144-146);
- the other mate was then aligned to the guide's taxon T (its own anchor of T, or looked for on T's gene and the gene
  the clade has past its end) and, on one gene or two neighbouring genes, completed the pair;
- the other mate's own best alignment, to a taxon U, was not weighed: the pair's Bitscore became the guide's plus the
  rescued mate's, and its MAPQ was taken against the next candidate, often U's lone mate, so it doubled.

Pairs whose mates stayed on two genes were written by `WriteSplitMates`, taking `ConsensusOfRead`'s taxon: a part
without a candidate of a taxon was no evidence against it (`ReadConsensus.h`, lines 88-100), and each mate got the
higher of its own and the consensus MAPQ.

## 2. The vote (`VoteMates`, `src/Core/MateGuidance.h`)

When the best candidate after joining and sorting has one mate only (one mate did not align, the mates prefer two
taxa, or they lie on two genes that are no neighbours), the two mates vote:
1. **Every taxon of the candidates** is scored by the sum of the two mates' best scores on it, `Score(2, 3, 1, 2)`
   as in a pair's Bitscore.
2. **A mate without an alignment of a taxon is looked for there** (`RescueMate`, guidance's own search, now a function
   of its own): its own anchor of the taxon, else the 12-mer diagonal on the other mate's gene where the fragment can
   reach and on the clade's neighbouring gene. Taxa are tried in the order of their known mate's score, at most 4
   (`kVoteMaxRescues`), and only while finding the mate (at most a perfect alignment) could still change the best or
   the second sum.
3. **A mate not found on a taxon but aligned elsewhere scores what an alignment of its read at the identity floor
   would** (`-a`, 0.9: 1.5 a base, 225 for 150 bases), at most its own best score. A mate aligned nowhere scores 0, as
   in a pair's Bitscore. Why: section 4.2; 0 lets a congener whose two mates align poorly beat the right taxon whose
   reference lacks one mate's part.
4. **The taxon of the highest sum wins**; of equal sums, the one of the better candidate (deterministic). Its **MAPQ**
   is `MAPQv2` of its sum against the next taxon's, over the mates aligned to either of the two (a mate aligned to
   neither adds the same to both). Mates that disagree get a low MAPQ; a mate found on the winner no longer doubles
   the pair against its competitors, since the competitors' mates are looked for too.
5. **Each mate is written at its best alignment of the winner**, its own or the one found, with the vote's MAPQ.
   A mate without an alignment of the winner, or whose own alignment of another taxon is **clearly better** (MAPQ
   20 or more between the two, `kVoteClearlyBetterMapq = kGuideMinMapq`: ~12 mismatches on 150 bases, or ~60
   bases clipped at a gene's end), keeps its own alignment with MAPQ 0, so that it counts for no taxon. A smaller
   difference is the strain's, and the pair, one organism, outweighs it.
6. A mate found on the winner completes the other mate's candidate where the two are one fragment (one gene in
   opposite orientations, or two neighbouring genes across their facing ends, as `JoinAlignmentPairs` pairs them);
   else it is a candidate of its own and the mates are written as split mates (`WriteSplitMates` with the vote's
   candidates and MAPQs, `PairVote` in `ReadConsensus.h`).

Unchanged: a pair whose best candidate already has both mates keeps today's ranking and MAPQ (against every other
candidate; the anchor recovery has already aligned each mate to the other's best genes). Pairs over two neighbouring
genes are paired by `JoinAlignmentPairs` as before; split pairs whose mates prefer different taxa get the vote's MAPQ.
Long reads keep `ConsensusOfRead`.

**The option.** `--mate_guidance consensus|guide|off` (`MateGuidanceMode.h`, `Options.h/.cpp`; `--no_mate_guidance`
= `off`); the options summary prints it, and the log has a line per sample: fragments voted on, with the mates' best
alignments on two taxa, moved to another taxon than the best candidate's, under MAPQ 4, mates looked for and found,
mates that kept another taxon's alignment.

**Speed work inside.** `mate_guidance::BestDiagonal` now sorts the k-mers of the shorter of read and window and
filters the other's through a 4,096-bit table before the binary search. Exact: a unit test compares it with the old
function on 300 random cases (repeats, ties, reads longer than the window); `guide` gives the same records as 0.8.4.

**Determinism.** The vote depends on the fragment's candidates alone (thread-local scratch, a fixed order); one
sample on 1 and 4 threads gave the same 7,872 records and the same profile ([check_final.sh](check_final.sh)).

## 3. Data and method

- **Database:** the reduced GTDB r226 database `~/c4test/db/protal_db/database.protal` (761 species, built by 0.8.4 at
  `5369c99`, pe gate 4 stride, forest models trained on 150 bp `guide` alignments).
- **Simulated test points** (the build's independent test set, `~/c4test/db/work/test/points`): rl150_p10000,
  rl150_p100000, rl150_p1000000 (10 samples each, 150 bp) and rl250_p10000, rl250_p100000 (250 bp). The rl150
  points' reads had been removed after the build; [resim.sh](resim.sh) made them again with the build's
  `simulate_metagenomes` and the commands in their `run_params.tsv`: manifests and abundance matrices identical to
  the originals, and `guide` gives the original SAMs' records exactly (check 3 below).
- **Two databases for them:**
  - `fin`: the finished database (every species of the test set in it), its own models;
  - `tdb`: the build's training database (the held-out species absent, as the build tested its models) with the
    finished database's pe model (`--model model_logs/trained_model.xml`).
- **Mocks:** the six paired-end mock samples in `~/Mocks/reads` (D6331 DNA, D6331 HKU and Auckland, Treichel even70,
  hCom2 in vitro, hCom1a cecum; 1M pairs each), on the finished database, scored by
  `scripts/mock_benchmark/score_mocks.py score`, with each sample's truth taxids as `PROFILE_TRUTH`.
- **Truth per fragment** (simulations): a read's name is `<contig>-<n>`; the contig's genome comes from the build's
  `genome_contigs.tsv.gz`, its species from `genomes_simulated.tsv`, a record's species from `genome2tiid.tsv`. A
  fragment "counts" for the species of its records with MAPQ 4 or more (the profiler's cutoff).
- **Scripts** (all here):
  [run_all.sh](run_all.sh) ([run_sims.sh](run_sims.sh), [run_mocks.sh](run_mocks.sh)) runs every sample per mode;
  [analyze.sh](analyze.sh) makes the tables with [score_sims.py](score_sims.py) (calls), [compare_records.py](compare_records.py)
  (records), [changed_fragments.py](changed_fragments.py) (fragments that moved), [feature_shift.py](feature_shift.py)
  (the model's 90 features and its calls) and [aligner_counts.py](aligner_counts.py);
  [count_instructions.sh](count_instructions.sh) counts instructions (callgrind); [check_final.sh](check_final.sh)
  checks the final binary. Results in [results/](results/): `guide`, `consensus_floor` (the rule kept, the final
  default), `consensus_zero`, `consensus_neutral` and the `_clip` variants (section 4.2); `fin` the finished
  database, `tdb` the training database.
- **Variants of the unaligned mate's score** (section 4.2) were run with intermediate builds of this branch: the first
  vote scored a mate not found 0 (`zero`); a second build chose by an environment variable between `floor` (the rule
  kept) and `neutral` (its best score elsewhere, as `ConsensusOfRead` treats a part without a candidate), and a third
  could also count the bases clipped at a gene's end at the floor rate (`_clip`). The final binary's default run gives
  the `floor` build's records exactly (97,689 fragments compared, check 1).

Commands (WSL):

```bash
D=docs/claude/2026-10-15-mate-consensus
bash $D/resim.sh ~/mc/sims rl150_p10000 rl150_p100000 rl150_p1000000
bash $D/run_all.sh ~/mc/bin/protal ~/mc/runs        # guide and consensus: fin, tdb, mocks
bash $D/analyze.sh ~/mc/runs ~/mc/results consensus
bash $D/check_final.sh ~/mc/bin/protal ~/mc/runs ~/mc/check
R=~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922
bash $D/count_instructions.sh ~/mc/bin/protal ~/c4test/db/protal_db $R/SRR23926922_1.fq.gz $R/SRR23926922_2.fq.gz ~/mc/cg 200000
```

## 4. Results

### 4.1 The calls

Simulated test points (50 samples each; "rows" are the `.truth_annotated` tables, as the build counts its test set;
"calls" the profiles against every simulated species the database has, so FN includes species without a row)
([results/sims_fin.summary.tsv](results/sims_fin.summary.tsv), [results/sims_tdb.summary.tsv](results/sims_tdb.summary.tsv)):

| database | mode | calls TP / FP / FN | F1 | rows TP / FP / FN | rows F1 |
|---|---|---|---|---|---|
| finished | guide | 1668 / 10 / 410 | 0.8882 | 1668 / 10 / 31 | 0.9879 |
| finished | consensus | 1669 / 11 / 409 | 0.8882 | 1669 / 11 / 30 | 0.9879 |
| training (held-out absent) | guide | 1087 / 26 / 991 | 0.6813 | 1087 / 26 / 60 | 0.9619 |
| training (held-out absent) | consensus | 1089 / 25 / 989 | 0.6823 | 1089 / 25 / 58 | 0.9633 |

- **Per point** the two modes differ by one species at most: finished database rl150_p100000 +1 TP, rl250_p100000
  +1 FP; training database rl150_p100000 +1 TP, rl150_p1000000 +1 TP and −1 FP.
- **The build's own independent test set** is the training database's three rl150 points scored by rows. `guide`
  gives exactly the build's figures (TP 678, FP 18, FN 39, F1 0.9597); `consensus` gives TP 680, FP 17, FN 37
  (F1 0.9618).
- **Mocks** ([results/mocks_samples_by_run.tsv](results/mocks_samples_by_run.tsv)): the same calls in all six
  samples (D6331 DNA 10/0/5, Treichel 62/0/1, hCom2 53/0/64, D6331 HKU and Auckland 12/0/3 each; hCom1a 54 called).
  Abundances move by less than 0.0003 of L1. The unknown share moves by −0.0017 (hCom2) to +0.0013 (HKU).

### 4.2 The records, and how a missing mate should score

Fragments with a record in either run; a fragment counts for the species of its mates with MAPQ 4 or more; "right"
is the species it was simulated from ([results/records_*.tsv](results/), [results/changed_*.txt](results/)):

| | finished db: guide | consensus | training db: guide | consensus | mocks: consensus |
|---|---|---|---|---|---|
| fragments with records | 615,026 | | 483,780 | | 247,303 |
| fragments whose records differ | | 6,457 (1.05%) | | 4,425 (0.91%) | 4,909 (1.99%) |
| of them a mate's taxon | | 5,658 | | 3,237 | 4,477 |
| crossing MAPQ 4 | | 4,246 | | 2,880 | 2,454 |
| mates written on two taxa | 215 | 227 | 684 | 667 | 366 → 291 |
| counted mates on the right species | 1,029,456 | 1,032,954 | 665,698 | 667,224 | |
| counted mates on a wrong species | 8,971 | 9,024 | 152,494 | 152,739 | |
| fragments counted for the right species | 573,032 | 572,775 (−257) | 370,912 | 370,838 (−74) | |
| fragments counted for a wrong species | 5,796 | 5,777 (−19) | 92,835 | 92,440 (−395) | |

- **More mates are placed** (+3,500 counted mates on the right species, finished database): the vote looks for a
  missing mate on every competing taxon, not only when a mate is sure. Of 208,942 mates looked for on the finished
  database 29,047 were found (guidance: 16,241 found over 111,649 fragments guided).
- **Fragments that lost their right species** (310 on the finished database) are mostly borderline: 255 went from
  MAPQ 4 to 3. A mate found as well on the runner-up taxon (often the same partial alignment at a gene end) joins
  both sums, and `MAPQv2` falls with the ratio of the sums, as it does for any pair on one gene.
- **Fragments that gained a wrong species**: 42, against 61 that lost one (finished database), and 74 against 469
  (training database, where the held-out species' reads can only go to relatives).

**How a mate not found should score.** Three rules were run (each against `guide`; right / wrong are the changes of
fragments counted for the right and a wrong species):

| rule for a mate not found on a taxon | finished db: right / wrong | training db: right / wrong | wrong gained / lost (finished) |
|---|---|---|---|
| `zero`: 0, as in a pair's Bitscore | −175 / **+68** | −30 / −20 | 109 / 41 |
| `floor`: an alignment at the identity floor, at most its best (kept) | −257 / −19 | −74 / −395 | 42 / 61 |
| `neutral`: its best score elsewhere (no evidence) | −333 / −49 | −126 / −616 | 56 / 105 |
| `zero` + clipped bases at the floor rate | −1,297 / −192 | −517 / −568 | 118 / 310 |
| `floor` + clipped bases at the floor rate | −1,307 / −307 | −527 / −959 | 17 / 324 |

- **`zero` moved fragments to the wrong species.** Its 109 new wrong ones ([results/changed_fin_consensus_zero.tsv](results/changed_fin_consensus_zero.tsv))
  are mostly a sure mate on its own species whose other mate the species' reference lacks (a gene, or a gene's
  neighbour, missing, or a stretch too divergent in the reference genome), while a congener has both:
  - *Zymomonas sp020216885*: mate 1 150M on it (MAPQ 100 under `guide`), mate 2 only on *Z. mobilis* with 16
    mismatches. `zero`: 300 + 0 against 270 + 220, so both mates went to *Z. mobilis* with MAPQ 44. `floor`:
    300 + 220 against 490, the right species keeps mate 1.
  - *Veillonella_A sp015058545* and *Micromonospora* sp.: the same pattern.
- **`neutral` is the most lenient**: fewest wrong fragments on the training database, but the most right ones lost.
  A mate's absence then costs nothing, so a partial alignment at a gene's end (e.g. 28 bases, 56) loses to a
  congener where the mate is simply absent.
- **Counting clipped bases** (to undo the gene-end bias, audit L2) removed the most wrong fragments, but it gave every
  partial mate a large equal score on both competing taxa, which diluted the MAPQ of the other mate's margin: ~1,300
  right fragments fell under MAPQ 4 (finished database), five times `floor`'s loss. Rejected.
- **`floor` kept**: it removes `zero`'s failure, it keeps the right species' wrong-to-right balance, and its premise
  is plain: a mate the database cannot place is worth at most the poorest alignment protal accepts.
- **In calls the rules hardly differ**: finished database 1668-1669 TP, 10-12 FP; training database 1087-1089 TP,
  25-26 FP; mocks identical.

### 4.3 Features under models trained on `guide`

[feature_shift.py](feature_shift.py) compared the 90 features of the pe model per taxon and the model's probability
([results/shift_*_consensus_floor.*](results/)):

| | finished db | training db | mocks |
|---|---|---|---|
| taxa with a row in both | 2,091 | 1,903 | 409 |
| probability, mean change / mean absolute change | +0.0007 / 0.0022 | +0.0010 / 0.0030 | +0.0022 / 0.0115 |
| calls that changed | 4 (2 TP gained, 1 TP lost, 1 FP gained) | 3 (2 TP gained, 1 FP lost) | 0 |

- **Every changed call is at the knob's edge**: probabilities 0.45-0.63.
- **The features that moved most:**
  - the sample's read-length summaries, which change for every taxon of a sample because more partial mates are
    placed: `sample_read_length_sd` +1.2 bp (of ~40) on the simulations and +1.5 bp (of ~26) on the mocks,
    `aligned_length_ratio` by at most 0.003;
  - `low_mapq_share` of absent taxa +0.01;
  - `mate_lost_share`, `congener_fit_share`, `ancestry_congener_share` and `variant_sites_per_kb` of a few taxa.

  No feature moves a present or absent taxon's mean by more than a few percent.

**What a retrained model would need.** Nothing beyond a database built at this commit: `build_gtdb_database.py` runs
protal with its default, so a new build trains its models on `consensus` alignments. It should record
`mate_guidance consensus` in `build_settings.tsv`. A retraining was not run here. With 0-4 calls changed among
~2,000 taxa per set, all of them at the knob, the old models are not the limiting factor.

### 4.4 What the aligner counted

Summed over each run's samples ([results/aligner_counts.tsv](results/aligner_counts.tsv)):

| | finished db | training db | mocks |
|---|---|---|---|
| `guide`: fragments guided / mates found | 111,649 / 16,241 | 108,513 / 15,343 | 22,413 / 3,140 |
| `consensus`: fragments voted on | 149,294 | 132,720 | 52,859 |
| of them, mates' best alignments on two taxa | 298 | 876 | 471 |
| moved to another taxon than the best candidate's | 93 | 105 | 339 |
| vote MAPQ under 4 | 14,098 | 8,527 | 18,123 |
| mates looked for / found | 208,942 / 29,047 | 167,351 / 22,251 | 110,220 / 15,537 |
| mates that kept another taxon's alignment (MAPQ 0) | 227 | 667 | 294 |

Real reads disagree more often than simulated ones: mates on two taxa are 0.9% of the mocks' voted fragments and
0.2% of the finished-database simulations'.

### 4.5 Speed

Instructions of the paired-end alignment loop (the OpenMP body of `RunPairedEnd`, callgrind, one thread), the first
200,000 pairs of D6331 HKU on the finished database ([count_instructions.sh](count_instructions.sh),
[results/instructions.txt](results/instructions.txt)):

| `--mate_guidance` | instructions | against `off` | against `guide` | mates looked for / found |
|---|---|---|---|---|
| `off` | 6,532,876,348 | | | |
| `guide` | 6,600,694,651 | +1.04% | | 1,131 / 196 |
| `consensus` | 6,775,903,343 | +3.72% | **+2.65%** | 4,360 / 688 |

- **Where the vote's cost goes:** `VoteMates` takes 240.7 M instructions inclusive (3.55% of the loop). Of that,
  looking for mates (`RescueMate`) is 236.8 M, and the 12-mer diagonals in it (`BestDiagonal`) 152.0 M.
- **How often it runs:** 2,201 of the 200,000 fragments were voted on, and 35 of them had their mates' best
  alignments on two taxa.
- **Where the cost could fall:** 84% of the mates looked for are not found. Scanning the window without storing its
  k-mers, or a bound on the missing mate's score tighter than a perfect alignment, would cut it; neither was needed
  for the few percent asked.

Wall-clock times on the shared machine were no use (load 5-22): the same runs differed by more than the whole loop
between repeats.

### 4.6 Tests and checks

- **Unit tests:** 592 pass (`ctest`), among them new ones:
  - `MateGuidance.BothMatesVote`: ties, a mate not found, both mates on both taxa, the clearly-better rule, a mate
    the right taxon's reference lacks, a mate aligned nowhere, no vote for a pair;
  - `MateGuidance.TheDiagonalIsTheSameWhicheverSideIsSorted`;
  - the vote across neighbouring genes in `AcrossGenes.MateGuidanceLooksPastTheGenesEnd`.
- **End-to-end tests:** 182 pass with the default (`consensus`) on a fresh mini database. With `--mate_guidance
  guide` passed by a wrapper, all pass but one, which compares the binary's own path in a message (the wrapper's).
- **[check_final.sh](check_final.sh)** ([results/check_final.txt](results/check_final.txt)):
  1. the final default run gives the evaluation's `floor` records (97,689 fragments, 0 differ);
  2. 1 and 4 threads give the same records and profile;
  3. `guide` gives the build's own 0.8.4 records on the training database (447,059 fragments, 0 differ in taxon or
     MAPQ).

## 5. Default and open points

**Default `consensus`.**
- The calls are neutral to slightly better (simulations ±1-2 species, mocks identical) with models trained the old way.
- Fewer fragments count for a wrong species, and mates that disagree get a MAPQ that says so instead of a doubled one.
- The alignment loop costs 2.65% more instructions than with `guide`.
- `guide` stays selectable for comparison and for databases whose users want their models' exact alignments.

**Open.**
- **`build_settings.tsv`** must record `mate_guidance` (another session's work), and runs should warn when theirs
  differs from the database's, as for the WFA2 traceback.
- **MAPQ dilution:** a mate found equally on both taxa lowers the vote's MAPQ, as for any pair on one gene (255
  right fragments from MAPQ 4 to 3 on the finished database). A MAPQ over the mates that tell the taxa apart would
  keep them, but it would differ from how pairs on one gene are scored, so it was left out.
- **Pairs aligned together are not voted on.** A pair on T whose mate 1 fits U better, while mate 2 has no candidate
  of U, keeps today's MAPQ against U's lone mate. Voting there would look for mates on most pairs and cost far more.
- **The gene-end bias** (audit L2, clipped bases count 0) remains, in the vote as in every pair; the attempt to
  correct it inside the vote cost more than it gained (section 4.2).
- **Not run:** a retrained model on `consensus` alignments (section 4.3), and r226 at full scale.
- **Found on the way** (not changed here): with a map whose `PROFILE` names a file in `profiles/`, the first sample
  profiled loses its `.truth_annotated` ("Writing the truth annotation failed"). `AnnotateWithTruth` (Profiler.h)
  writes it before the profile has created the folder. [run_sims.sh](run_sims.sh) makes the folder first.

**Website.** The options page should list `--mate_guidance consensus|guide|off` (default `consensus`) in place of
`--no_mate_guidance` (kept as an alias of `off`), and its SAM section the vote's MAPQ for split pairs.

**Data.** WSL `~/mc/runs/{fin,tdb,mocks}/{guide,consensus_floor}` keep the SAMs, profiles and truth tables of the
two modes compared (`consensus_floor` is the final default's run), and `~/mc/results2` all the tables, the variants'
included. Removed: the other variants' runs, the rl150 reads (resim.sh remakes them) and the build trees.

## 6. After the review of the merged branch (2026-10-15)

The integration branch's code review found that the vote's MAPQ compared the winning taxon with the next taxon only:
a mate that fits two places of the winner equally well (a paralog, another gene) got a high MAPQ where `MAPQv1` had
given it 0. Since `126f800` on `f1-integration` each mate's MAPQ is also at most its margin within the winner (its best
alignment there against its next one at another place; `MateGuidance.h`, unit test in `test_SamRoundTrip.cpp`).
On the six paired-end mocks (`run_mocks.sh` with `consensus`, the binaries of `848c6d1` and `126f800`;
`compare_records.py`, then the calls) 13 of 247,303 fragments changed, their MAPQ only; none crossed the profiler's
cutoff of 4, and no taxon's probability moved. The numbers of section 4 stand.

The mode is also recorded with the database since the merge with the train/run guards: `--build` writes
`mate_guidance` (the mode given, else `consensus`) into `build_settings.tsv`, and runs follow it unless
`--mate_guidance` is given; a database that records none (built before 2026-10-15) keeps `guide`, which its models were
trained with ([running.md](../../running.md#what-a-run-takes-from-its-database-and-models)).
