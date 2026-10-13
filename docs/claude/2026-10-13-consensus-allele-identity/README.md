# *E. coli* K-12 in HiFi mocks: identity against the best strain allele, and a consensus tie-break left out

2026-10-13, branch `audit-fixes`. A follow-up of [2026-10-13-ecoli-hifi](../2026-10-13-ecoli-hifi/README.md), which found
why HiFi reads of *E. coli* K-12 are called *Escherichia sp004211955*. The user asked for two of its proposals:
- **Minimax consensus tie-break:** when every taxon of a long read loses to one, the taxon whose largest defeat is
  smallest takes the read.
- **Identity features against the best strain allele.**

**What was committed.**
- **Committed:** the identity features, as the new feature group `strain_identity`
  ([features.md](../../features.md#the-identities-against-the-best-strain-allele-strain_identity-2026-10-13)).
- **Not committed:** the tie-break. It was implemented and tested, and moved more reads away from the sample's species
  than towards them (section 3); its patch is [minimax_tiebreak.patch](minimax_tiebreak.patch).

**Code.**
- **Tested:** both changes on `49cd583`.
- **Committed:** on top of `19ccc0e`, which adds the opt-in gene grid and changes no default output.

**Data.** The reduced r226 database of the mock local test ([2026-10-13-mock-local-test](../2026-10-13-mock-local-test/README.md):
761 species), built again from fresh inputs, and the test's twelve real samples in `~/Mocks`.

## 1. What changed

**Identity against the best allele** (`Profiler.h`).
- **Per record:** the profiler already scores each kept record against its copy's strain alleles (`NoteAlleles`, for
  `allele_explained_share` and `allele_identity_gain`). It now also returns that record's gain: the differences its
  best allele explains, less the allele's own edits the read lacks, per aligned column.
- **Per read:** the gain is carried with the record to its gene (`SamAddition`, `ReadEvidence`). There each read gets a
  second identity, its identity plus the gain, at most 1.
- **Two new features:**
  - `allele_identity`: `identity` computed from those identities;
  - `allele_top_identity`: `top_identity` computed from them.
- **Without alleles:** for a species without alleles, or a database without the table, they equal `identity` and
  `top_identity`.
- **Unchanged:** the depth and the strain MSAs keep the identities against the representative.
- **Feature group:** `strain_identity` (`scripts/model_features.py`), in the default set: 93 features.
  - **Older databases:** their models do not name the two features, so they score as before.
  - **Older training tables:** these lack the two columns. They train with the set without the group, which
    `--features auto` still offers.

**The tie-break** (`ReadConsensus.h`, not committed).
- **The current rule:** taxa are compared on the parts (genes) both have candidates on, by the sum of their bitscores.
  The consensus is a taxon nobody beats. When every taxon is beaten, it is the one with candidates on the most parts.
- **The patch:** in that case it takes the taxon whose largest defeat is smallest, among those with candidates on
  enough parts.
- **Not the default:** ReadConsensus.h's comment records that it was tried.

## 2. The comparison

**Four reduced databases** were built from the same inputs, seed and settings (as the mock local test's, scaled to 4
cores; [build_db.sh](build_db.sh), [run_all.sh](run_all.sh)):

| Build | protal | Models trained on |
|---|---|---|
| `base` | `49cd583` | the default set of `49cd583` (91 features) |
| `identity` | + the allele identities | the default set with `strain_identity` (93) |
| `consensus` | + both changes | the default set without `strain_identity` (91): the tie-break alone |
| `both` | + both changes | the default set with `strain_identity` (93) |

Each build took 10-22 min on the shared machine. The mock samples were profiled with each build's protal and database
([run_mocks.sh](run_mocks.sh)), and compared with [compare.py](compare.py)
([results/compare_runs.txt](results/compare_runs.txt)).

**The rebuild is exact.** `base` gives the same calls and the same probabilities, to the digit, as the original run
`0.8.3_4dfc549` of the mock test. The inputs were made again (`make_inputs.sh`; 1,028 of 1,038 genomes again), so the
differences below are the changes'.

### The identity features

**Mock samples**, where anything changed (all other samples: the same calls):

| Sample | `base`: TP/FP/FN, F1 | `identity`: TP/FP/FN, F1 | Change |
|---|---|---|---|
| ATCC MSA-1003 HiFi | 17/1/3, 0.894 | **19/0/1, 0.974** | *E. coli* called (P 0.12 → 0.87), *Helicobacter pylori* called (0.3% of cells), *sp004211955* no longer called (0.55 → 0.35) |
| Zymo D6331 HiFi (WGA) | 13/2/2, 0.866 | **14/1/1, 0.933** | *E. coli* called (0.12 → 0.97), *sp004211955* no longer called (0.63 → 0.49) |
| D6331 ONT | 14/0/1, 0.965 | 14/0/1, 0.965 | *E. coli* 0.59 → 0.63 |
| D6331 pe, HKU / Auckland | 12/0/3 each | the same | *E. coli* 0.80 → 0.92, 0.82 → 0.94 |
| hCom2-like in vitro (pe) | 52 of 117 | 53 of 117 | *Phocaeicola vulgatus* called |

Abundances are unchanged: L1 0.462 → 0.458 ATCC, and in D6331 HiFi 0.765 → 0.624 only because *E. coli* is now
reported.

**The features themselves** (`--profile_truth` of the `identity` run, [dump_features.sh](dump_features.sh),
[results/feature_comparison_identity.tsv](results/feature_comparison_identity.tsv)):

| Sample | Taxon | P | `top_identity` | `allele_top_identity` | `identity` | `allele_identity` | `allele_explained_share` |
|---|---|---|---|---|---|---|---|
| ATCC HiFi | *E. coli* | 0.87 | 0.9957 | **1.0** | 0.980 | 0.996 | 0.84 |
| ATCC HiFi | *sp004211955* | 0.35 | 0.9984 | 0.9984 | 0.987 | 0.987 | 0 |
| ATCC HiFi | *H. pylori* | 0.94 | 0.981 | **1.0** | 0.968 | 0.991 | 0.79 |
| ATCC HiFi | true positives (median of 19) | 1.0 | 1.0 | 1.0 | 0.993 | 0.997 | 0.33 |
| D6331 HiFi | *E. coli* | 0.97 | 0.9957 | **1.0** | 0.979 | 0.995 | 0.80 |
| D6331 ONT | *E. coli* | 0.63 | 0.994 | 1.0 | 0.968 | 0.982 | 0.54 |

- ***E. coli*:** K-12's best reads match a known *E. coli* genome exactly, so against the best allele it looks like
  every other present species.
- ***sp004211955*:** the congener has no alleles, so its identities stay as they were.

**The build's own test sets** (simulated; `model_logs/summary.txt`, F1 at the build's knob):

| Read type | `base` | `identity` | `both` |
|---|---|---|---|
| pe, species held out | 0.9736 | 0.9748 | 0.9748 |
| pe, independent test set | 0.9631 | 0.9618 | 0.9618 |
| HiFi, species held out | 0.9734 | 0.9788 | 0.9781 |
| HiFi, independent test set | 0.9680 | 0.9680 | 0.9680 |
| ONT, species held out | 0.9701 | 0.9687 | 0.9713 |
| ONT, independent test set | 0.9620 | 0.9602 | 0.9583 |

All within noise:
- **pe test set:** 2 more false positives in 1,304 taxa.
- **ONT test set:** 1 more in 321.
- **HiFi, species held out:** 5 fewer false negatives and 3 fewer false positives in 903.

The simulated strains are drawn from genomes that are never alleles, so the training sees few strains as close to
an allele as K-12 is; the gain shows on real strains.

### The consensus tie-break

**At the read level** the change touches few reads, and moves more of them away from the sample's species than
towards it.
- **Method:** [consensus_changes.py](consensus_changes.py) takes each read's consensus taxon (the taxon most of its
  records name, `ZR:i:2` left out) in two SAMs of the same sample, and sorts the reads that changed by whether they
  moved to or from a species of the truth.
- **Variants:** a scratch build of the same protal had the rule chosen at run time
  ([ReadConsensus_variants.h](ReadConsensus_variants.h)); `old` reproduces `base`'s SAMs exactly.
- **Rules tried:** besides the requested tie-break, a per-gene vote: taxa compared by the common parts each wins, the
  sum breaking a tie, which the example read below suggested.

| Rule | ATCC HiFi: to truth / away | D6331 HiFi | D6331 ONT |
|---|---|---|---|
| minimax tie-break (the request) | 2 / 8 | 0 / 55 | 2 / 19 |
| per-gene vote | 2 / 16 | 1 / 10 | 4 / 7 |
| per-gene vote with the tie-break | 3 / 29 | 1 / 92 | 10 / 48 |

Of ~20,000-25,000 reads per sample. The main moves under the tie-break:
- *Veillonella rogosae* → *V. sp900549845*: 40 reads in D6331 HiFi, 16 in ONT;
- *E. coli* → *E. fergusonii* or *sp004211955*: 8 in ATCC HiFi, 10 in D6331 HiFi.

**Why the tie-break favours relatives.** A relative with candidates on few of the read's parts is compared with the
others on those parts only, so its defeats are small.
- **The example:** read SRR13128012.101382 of D6331 (*V. rogosae*) has two genes.
  - *V. sp900549845* has a candidate on one of them, where it beats *rogosae* by 45 points (9 mismatches).
  - *rogosae* beats every taxon with both genes, by 190 points or more.
  - In the resulting cycle, *sp900549845*'s largest defeat (5 points, against *V. dispar*) is the smallest, so it
    takes the read. The current rule gives it to *rogosae*, which has candidates on both.

**Why the K-12 reads are not saved by it.** Read SRR11606871.13443 of the earlier report, its segments' candidates
printed by a debug build ([results/read13443_segments.tsv](results/read13443_segments.tsv)):
- **Same footing:** all taxa's alignments cover the same read bases, so no copy is favoured by its length.
- **Gene 51 decides:** *E. coli* has 34 mismatches there, *E. fergusonii* 12. The earlier report could not see this:
  `ZA` lists only candidates within 5 edits of the best.
- **The cycle:**
  - *E. coli* loses to *fergusonii* on their 7 common genes by 25 points, although it wins 3 genes to 2;
  - *fergusonii* loses to *E. whittamii* by 25 points; *whittamii* has no candidate on gene 51, so that gene counts
    for nothing against it;
  - *whittamii* loses to *E. coli* by 65 points.
- **The outcome:** *E. coli* and *fergusonii* tie for the smallest largest defeat, and *fergusonii*, with candidates on
  more genes, takes the read.
- **Before:** *sp004211955*, which had candidates on the most genes. The consensus also steers which genes are looked
  for next to the read's genes, so the two rules end with somewhat different candidates: the old one after two
  passes, the tie-break after three.

These reads are not an artifact of the fallback: one gene of K-12 that is nearer a congener outweighs the rest of the
read. The per-gene vote would give this read to *E. coli*, but on the samples as a whole it too moves more reads away
from the truth than towards it.

**At the profile level** the tie-break adds nothing to the identity features:
- **The `consensus` build** (the tie-break alone):
  - **HiFi:** the wrong *sp004211955* call is gone, as its reads go to *E. fergusonii* instead. ATCC 17/0/3, F1 0.918;
    D6331 13/1/2, 0.896. *E. coli* is still missed (P 0.12).
  - **ONT:** a trace of *Klebsiella sp004342285* (6 × 10⁻⁶) is called: 14/1/1, F1 0.933.
  - **Its test sets:** HiFi 0.9680 → 0.9660, the others as `base`.
- **`both` against `identity`:**
  - **D6331 HiFi:** *sp004211955* is still called (P 0.53).
  - **ONT:** the same trace of *Klebsiella*.
  - **F1:** 0.903 and 0.933, against `identity`'s 0.933 and 0.965.

**A correction to the earlier report.** Its `lost_reads.py` called a read "cycle" when *E. coli*'s `ZA` edits were at
most those of the record written on most of its genes. Genes where *E. coli* is more than 5 edits worse are missing
from `ZA`, so reads such as 13443 were counted as cycles that the consensus decided on a gene invisible there. Its
"locus" and "mixed" classes are not affected. Its loss of 18% of *E. coli*'s reads to congeners is unchanged by every
rule here.

## 3. Conclusions

- **The identity features fix the calls.** Both HiFi misses are called, and both wrong congener calls are gone, with
  the simulated test sets within noise.
- **They need a rebuild.** A database's models must be trained with the new group. The next r226 build takes them by
  default. To check there: the ATCC MSA-1003 and D6331 HiFi runs (`download_mocks.py --sets
  atcc1003_hifi,d6331_hifi`).
- **The consensus stays as it is.** Neither the tie-break nor a per-gene vote improves the read assignment on these
  samples.
- **The remaining loss is biology.** K-12 loses a fifth of its long reads to congeners, through genes nearer a
  congener's copies (the earlier report's "locus" reads and reads like 13443). With the allele identities this no
  longer costs the call.

## Tests

**`new` tree** (`49cd583` + both changes + the map-folder test fix below):
- **C++ unit tests:** all pass, including:
  - `StrainAlleles.TheIdentityAgainstTheBestAllele`, new: a read with its allele's two differences reaches 1.0, one
    with another difference keeps its identity, and a species without alleles gets the plain identities;
  - `ReadConsensus.ACycleGoesToTheTaxonWithTheSmallestDefeat`, part of the left-out patch.
- **Other suites:** e2e 177 tests, the model tests 45, the mini-db and script tests 126.

**The map-folder test.** `SampleMap.ARelativeOutputFolderIsUsedOnce` had failed since `d865681` (no empty output
folders): it still expected reading a map to make `profiles/`. Its fix, found uncommitted in the checkout, is committed
with this change.

**The commit** was checked from a `git archive` of it: see the commit message.

## Files

| File | What it does |
|---|---|
| [build_db.sh](build_db.sh) | the reduced database, by one protal clone, with more build options |
| [run_mocks.sh](run_mocks.sh) | protal of one clone on the twelve mock samples, and the scores |
| [run_all.sh](run_all.sh) | the four builds and their mock runs |
| [compare.py](compare.py) | runs side by side: scores, *Escherichia* probabilities, *E. coli*'s lost long reads |
| [consensus_changes.py](consensus_changes.py) | the reads whose consensus taxon moved between two SAMs |
| [dump_features.sh](dump_features.sh) | a run's long-read samples re-profiled with `--profile_truth` |
| [ReadConsensus_variants.h](ReadConsensus_variants.h) | the consensus rules tried, chosen by `PROTAL_CONSENSUS` (scratch, not in the build) |
| [minimax_tiebreak.patch](minimax_tiebreak.patch) | the tie-break and its unit test, as implemented on `49cd583` |
| [results/](results/) | the comparison, the feature table, read 13443's candidates |

**Kept in WSL:**
- `~/idtest`: the clones, inputs and four builds, 9 GB. Delete it when done.
- `~/Mocks/idc_*`: the four runs.
