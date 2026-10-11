# Strain false negatives: strain alleles, indexed alleles and more genomes per species

**Date:** 2026-10-11. **Code:** protal `8d8c89b` (WSL `~/b23`, Release; sources checked against HEAD), branch
`audit-fixes`. **Machine:** WSL, 4 cores (`taskset -c 0-3 nice -n 5`), shared with another session's compile.

**Data.**
- **The world** ([select_strain_world.py](select_strain_world.py), [selection.txt](selection.txt)): GTDB r226,
  bacteria. 30 focal species with 24 or more good genomes besides the representative (CheckM2 >= 90% complete,
  <= 5% contaminated), 10 per stratum of GTDB's mean intra-species ANI (tight >= 99.0%, mid 97.5-99.0%,
  wide < 97.5%), each from its own genus of 2-30 species. Every other species of those genera comes in as its
  representative. That makes 449 species and 1,520 genomes.
  - Each focal species has up to 40 other genomes, dealt by hash into an **allele half** (may enter a database's full
    reference) and a **test half** (never in any database; reads are simulated from it) ([halves.tsv](halves.tsv)).
  - The marker genes were streamed out of the packed `bac120_marker_genes_all_r226.tar.gz` in WSL `~/GTDB/r226` by
    the real-ancestry report's `extract_subset.py` ([run_extract.sh](run_extract.sh); 14 min). Everything lives in
    WSL `~/strainfn`.
- **The test strains:** 532, classed by their marker genes' 21-mer containment ANI to their representative
  ([evaluate.py](evaluate.py), [divergence.tsv](divergence.tsv)):

  | class | strains | median ANI to the rep | to the nearest of 5 allele-half genomes | of 20 |
  |---|---|---|---|---|
  | easy (>= 99.5%) | 198 | 99.75% | 99.94% | 99.98% |
  | medium (98.5-99.5%) | 248 | 98.96% | 99.32% | 99.47% |
  | hard (< 98.5%) | 86 | 98.27% | 99.43% | 99.90% |

- **The reads** ([simulate_strain_reads.py](simulate_strain_reads.py), [truth.tsv](truth.tsv)): paired 150 bp reads
  drawn from the marker genes only (insert N(350, 50), substitutions 0.1-0.6% along the read). There are 20 samples
  per set. Sample r holds the test strain of rank r of every focal species, plus 15 congener representatives at
  20-200 pairs.
  - **deep:** every strain at 60 pairs (read level).
  - **shallow:** each strain at 1, 2, 3, 4, 6, 8, 12, 16, 24 or 32 pairs, cycled (calls). That gives 600 strain
    instances per variant.
- **Databases** ([make_variants.sh](make_variants.sh), [build_lines.txt](build_lines.txt)): one conversion
  (`gtdb_to_protal_db.py`), then one `protal --build` per variant. All variants share the same representatives and
  full reference (the test half removed). They differ only in:

  | variant | build options | allele genomes per focal species |
  |---|---|---|
  | `reps` | `--strain_alleles 0` | (20, unused) |
  | `a4` (the default) | `--strain_alleles 4` | 20 |
  | `a4i01` | `--index_alleles 0.01` (as r226 v22/v23) | 20 |
  | `a4i0` | `--index_alleles 0` (every allele indexed) | 20 |
  | `a16` | `--strain_alleles 16` | 20 |
  | `a16i0` | `--strain_alleles 16 --index_alleles 0` | 20 |
  | `a16i0n5` | `--strain_alleles 16 --index_alleles 0` | 5 |
  | `oracle` | `--strain_alleles 16 --index_alleles 0` | all 40, the test strains too ([run_oracle.sh](run_oracle.sh)) |

- **Runs** ([run_samples.sh](run_samples.sh), [run_all.sh](run_all.sh)): one protal run per variant over all 40
  samples (`--no_strains --no_qcmsa`). Every run uses the same presence model: the paired-end GBM of the small real
  world, seed 1 (`~/realworld_sp_s1`, `371ac84`, trained with 4 alleles from half the genomes). Run-time variants on
  `a4`: `--no_allele_scores` (`a4ns`), `-c 6` (`a4c6`), `-w 16 -s 256` (`a4w16`); on `a16i0`: `-c 6`.
- **Evaluation:** [evaluate.py](evaluate.py) ([eval_all.txt](eval_all.txt), [eval_oracle.txt](eval_oracle.txt),
  `eval_*_reads.tsv`, `eval_*_calls.tsv`), [species_table.py](species_table.py) ([species_table.txt](species_table.txt)),
  [flips.py](flips.py) ([flips.txt](flips.txt): paired calls, sign test), [collect_costs.sh](collect_costs.sh)
  ([costs.tsv](costs.tsv), [run_costs.tsv](run_costs.tsv)).

**Questions (the user's).**
1. Do the extra alleles per species (intra-species variants) lead to more alignments and detections, in easy and
   hard cases?
2. Are the extra keys in the k-mer index (`--index_alleles`) needed for this?
3. Would more genome representatives per species at build time help, and what would they cost in memory?
4. Overall: how can the false-negative rate for strains of a known species be reduced?

## Summary

- **The alleles do not add alignments; they put alignments on the right species and above MAPQ 4.** Even with the
  representatives alone, 98.4-99.9% of a strain's pairs align somewhere. For hard strains, the alleles move the share
  counted for their species (right species, MAPQ >= 4) from 0.677 to 0.788, and with 16 indexed alleles to 0.810.
  Medium strains go from 0.715 to 0.764-0.774. Easy strains stay at 0.921-0.928.
- **They do add detections, mostly of hard strains.** Called test strains, over all depths:

  | | all (600) | easy | medium | hard | 1-2 pairs | FP (20 samples) |
  |---|---|---|---|---|---|---|
  | `reps` | 0.830 | 0.973 | 0.796 | 0.604 | 0.658 | 5 |
  | `a4` (default) | 0.892 | 0.982 | 0.844 | 0.822 | 0.758 | 2 |
  | `a16i0` | 0.903 | 0.982 | 0.851 | 0.871 | 0.792 | 3 |
  | `oracle` | 0.910 | 0.982 | 0.862 | 0.881 | 0.800 | 1 |

  `reps` to `a4`: 37 strains gained, 0 lost (p = 1e-11). Most of that comes from the alleles' **features**
  (`--no_allele_scores` on `a4`: 30 gained, 0 lost). The **scores** add 9 gained, 2 lost (p = 0.07), cut the pairs
  counted on another species by a fifth to a quarter, and bring FP from 5 to 2.
- **The extra index keys are not needed.** `a4` to `a4i01`: 6 gained, 2 lost (p = 0.29). `a16` to `a16i0`: 4 gained,
  2 lost (p = 0.69). Hard strains' counted pairs rise by about one point. In return the index grows by 2.5-5.8%
  here, and at r226 v22/v23 by 115M entries (database +1.1 GB, build memory after pass 2 +2.8 GB). Seeding is not
  where strains are lost.
- **More genomes per species help only a little beyond the default.** 16 alleles instead of 4 from the same 20
  genomes: 6 gained, 1 lost. 20 allele genomes instead of 5: 3 gained, 0 lost, and FP 6 to 3. Even the oracle,
  which includes the test strains' own genomes, gains only 4 more strains than `a16i0` (5 gained, 1 lost, p = 0.22).
  **The ceiling of more representatives is about +0.007 here.** Adding the genomes as references of their own
  would cost about +155% index entries at r226 (16 per species; section 6) and would split MAPQ between copies of
  the same species. As alleles, all of GTDB's genomes with 16 alleles per gene would cost about +0.9 GB of memory
  without indexing, or about +2.5-3 GB with it.
- **What is left is evidence-limited or congener-ambiguous, not missing references.** Of the 58 strains `a16i0`
  misses:
  - 21 have their right-species pairs only at MAPQ < 4;
  - 19 have a single counted pair;
  - 15 have 2+ counted pairs that the model rejects;
  - 3 aligned elsewhere or not at all.

  Two species account for 30 of the 58: *Dorea_A amylophila* (6/20 called) and *Priestia megaterium_D* (4/20). There,
  91-94% of right-species pairs fall below MAPQ 4 (`reps`) because a congener's marker genes are nearly identical;
  the oracle calls 7/20 and 5/20.
- **Other options:**
  - `-c 6` and `-w 16 -s 256` change nothing (+1 / −1 strain).
  - Knob 0.35 instead of 0.5 calls +1.8% (`a4`: 0.892 to 0.910) for FP 2 to 5 per 20 samples. This world's FP are
    too easy to judge that trade (section 5).
  - Counting MAPQ < 4 pairs would give 29 more missed strains 2+ pairs (17 to 46, `a4`), but would also give 2+
    pairs to 14× more absent species (10 to 141): not as a count.

## 1. Where a strain's reads go (deep samples, 60 pairs per strain)

Share of the strains' pairs ([eval_all.txt](eval_all.txt)); the best record of each pair:

| variant | class | aligned | right species | counted (MAPQ >= 4) | right, MAPQ < 4 | counted on another species | identity |
|---|---|---|---|---|---|---|---|
| `reps` | easy | 0.996 | 0.978 | 0.921 | 0.057 | 0.0044 | 0.9943 |
| `reps` | medium | 0.999 | 0.859 | 0.715 | 0.144 | 0.0258 | 0.9867 |
| `reps` | hard | 0.984 | 0.830 | 0.677 | 0.153 | 0.0513 | 0.9766 |
| `a4` | medium | 0.999 | 0.916 | 0.764 | 0.152 | 0.0208 | 0.9863 |
| `a4` | hard | 0.984 | 0.893 | 0.788 | 0.105 | 0.0390 | 0.9760 |
| `a4i01` | hard | 0.984 | 0.904 | 0.799 | 0.106 | 0.0318 | 0.9756 |
| `a16` | hard | 0.984 | 0.897 | 0.798 | 0.099 | 0.0377 | 0.9760 |
| `a16i0` | medium | 0.999 | 0.921 | 0.774 | 0.147 | 0.0187 | 0.9862 |
| `a16i0` | hard | 0.984 | 0.908 | 0.810 | 0.099 | 0.0298 | 0.9755 |
| `oracle` | medium | 0.999 | 0.925 | 0.781 | 0.144 | 0.0184 | 0.9862 |
| `oracle` | hard | 0.984 | 0.915 | 0.819 | 0.096 | 0.0273 | 0.9754 |

- **Aligned does not move** (identical in every variant). The representative's index already seeds nearly every
  pair of a hard strain at 98.3% ANI. The allele k-mers in the index (`a4i01`, `a16i0`) take a few pairs from
  congeners (counted on another species 0.039 to 0.030-0.032 for hard strains), not from "unaligned".
- **The scores move pairs to the right species.** Hard strains: 0.830 to 0.893-0.915 right, and counted on another
  species 0.051 to 0.027-0.039. `--no_allele_scores` gives `reps`' numbers exactly, since the alignments are the
  same.
- **About 15% of medium strains' pairs land on the right species under MAPQ 4 in every variant.** These pairs sit
  in genes that a congener shares nearly base for base. Mostly they come from two species (section 4); for them no
  allele can help.

## 2. Calls (shallow samples)

Share of test strains whose species is called (knob 0.5), by the strain's pairs ([eval_all.txt](eval_all.txt)):

| variant | class | 1-2 pairs | 3-6 | 8-16 | 24-32 | all | FP (20 samples) |
|---|---|---|---|---|---|---|---|
| `reps` | easy | 0.880 | 1.000 | 1.000 | 1.000 | 0.973 | |
| `reps` | medium | 0.577 | 0.780 | 0.871 | 0.911 | 0.796 | |
| `reps` | hard | 0.278 | 0.588 | 0.731 | 0.739 | 0.604 | |
| `reps` | all | 0.658 | 0.822 | 0.900 | 0.908 | 0.830 | 5 |
| `a4` | easy | 0.920 | 1.000 | 1.000 | 1.000 | 0.982 | |
| `a4` | medium | 0.673 | 0.829 | 0.894 | 0.946 | 0.844 | |
| `a4` | hard | 0.556 | 0.824 | 0.885 | 0.957 | 0.822 | |
| `a4` | all | 0.758 | 0.889 | 0.933 | 0.967 | 0.892 | 2 |
| `a16i0` | hard | 0.611 | 0.853 | 0.962 | 1.000 | 0.871 | |
| `a16i0` | all | 0.792 | 0.894 | 0.944 | 0.967 | 0.903 | 3 |
| `oracle` | hard | 0.611 | 0.882 | 0.962 | 1.000 | 0.881 | |
| `oracle` | all | 0.800 | 0.911 | 0.944 | 0.967 | 0.910 | 1 |

All variants, over all depths ([flips.txt](flips.txt) for the pairs):

| variant | all | hard | FP | vs the row it is compared with: gained / lost (sign test) |
|---|---|---|---|---|
| `reps` | 0.830 | 0.604 | 5 | |
| `a4ns` (`a4` without scores) | 0.880 | 0.812 | 5 | vs `reps`: 30 / 0 (2e-9) |
| `a4` | 0.892 | 0.822 | 2 | vs `a4ns`: 9 / 2 (0.07); vs `reps`: 37 / 0 (1e-11) |
| `a4i01` | 0.898 | 0.851 | 3 | vs `a4`: 6 / 2 (0.29) |
| `a4i0` | 0.897 | 0.851 | 3 | vs `a4`: 6 / 3 (0.51) |
| `a16` | 0.900 | 0.842 | 2 | vs `a4`: 6 / 1 (0.13) |
| `a16i0` | 0.903 | 0.871 | 3 | vs `a16`: 4 / 2 (0.69) |
| `a16i0n5` | 0.898 | 0.861 | 6 | vs `a16i0`: 0 / 3 (0.25) |
| `oracle` | 0.910 | 0.881 | 1 | vs `a16i0`: 5 / 1 (0.22); vs `a4`: 13 / 2 (0.007) |
| `a4c6` | 0.893 | 0.832 | 3 | vs `a4`: 3 / 2 (1) |
| `a4w16` | 0.890 | 0.822 | 2 | vs `a4`: 0 / 1 (1) |
| `a16i0c6` | 0.905 | 0.871 | 4 | vs `a16i0`: 1 / 0 (1) |

- **Easy strains** have nothing to gain: 0.973 to 0.982. Their 4 remaining misses with `a4` are all strains with 1-2
  pairs.
- **Hard strains** gain the most, and the gain grows with depth: at 24-32 pairs 0.739 to 0.957-1.000. With the
  representatives alone a hard strain at 24-32 pairs is missed one time in four.
- **The features carry most of the call gain.** The model was trained with alleles. Without them (`reps`), a species
  whose reads are 98% identical and unexplained looks like a novel congener's. With them the same reads are a known
  strain's.
- **Where an allele-half genome is nearer than the representative** (by 0.3% ANI or more; 223 of 532 strains), hard
  strains go from 0.631 to 0.845 (`a4`) and 0.881 (`a16i0`; 84 instances). Hard strains with no nearer genome go from
  0.471 to 0.706 and 0.824 (17 instances only).

## 3. Why the remaining strains are missed

The fate of each missed strain's pairs (shallow samples):

| variant | missed | none aligned | aligned to other species only | right species, all MAPQ < 4 | 1 pair counted | 2+ pairs counted (model said no) |
|---|---|---|---|---|---|---|
| `reps` | 102 | 1 | 3 | 41 | 22 | 35 |
| `a4` | 65 | 1 | 1 | 24 | 22 | 17 |
| `a16i0` | 58 | 1 | 2 | 21 | 19 | 15 |
| `oracle` | 54 | 1 | 1 | 21 | 18 | 13 |

The alleles roughly halve two groups: strains whose pairs all fall under MAPQ 4 (41 to 21) and strains the model
rejects with 2+ counted pairs (35 to 15). The group with a single counted pair does not shrink (22 to 19). A genome
nearer than the allele half's (the oracle) changes each group by one or two.

## 4. By species

[species_table.txt](species_table.txt) (called of 20, shallow samples; right-species pairs under MAPQ 4 from `reps`'
deep samples):

| species | stratum | median ANI to rep | right pairs under MAPQ 4 | `reps` | `a4` | `a16i0` | `oracle` |
|---|---|---|---|---|---|---|---|
| *Dorea_A amylophila* | mid | 0.9925 | 0.938 | 5 | 5 | 6 | 7 |
| *Priestia megaterium_D* | wide | 0.9895 | 0.906 | 1 | 5 | 4 | 5 |
| *Yersinia mollaretii* | wide | 0.9813 | 0.198 | 8 | 16 | 16 | 16 |
| *Lautropia mirabilis* | wide | 0.9801 | 0.371 | 9 | 14 | 16 | 17 |
| *Elizabethkingia miricola* | wide | 0.9860 | 0.357 | 14 | 17 | 17 | 17 |
| *Finegoldia magna_H* | wide | 0.9878 | 0.093 | 14 | 16 | 18 | 19 |
| *Metapseudomonas otitidis* | wide | 0.9875 | 0.017 | 15 | 18 | 20 | 20 |
| *Actinobacillus equuli* | wide | 0.9818 | 0.344 | 16 | 17 | 17 | 17 |
| *Phil1 sp001940855* | wide | 0.9852 | 0.152 | 16 | 19 | 19 | 19 |
| *Weissella cibaria* | mid | 0.9944 | 0.212 | 17 | 20 | 20 | 20 |
| *Proteus faecis* | tight | 0.9978 | 0.222 | 18 | 20 | 20 | 20 |
| the other 19 species | | | | 365 | 368 | 369 | 369 |

- **Two species hold half of `a16i0`'s misses.** In both, over 90% of right-species pairs fall below MAPQ 4: a
  congener's marker genes match them about as well as their own representative's (GTDB's *Priestia megaterium* and
  *megaterium_D*, two *Dorea_A* species). No reference sequence of the species can fix that: the oracle gains one
  strain each. Only the reads from the few discriminating sites can (longer reads, more depth).
- **The wide species gain most from alleles:** *Yersinia mollaretii* 8 to 16, *Lautropia mirabilis* 9 to 16-17.

## 5. Other options tried

- **More candidates (`-c 6`) and deeper seeding (`-w 16 -s 256`)** change no call beyond ±1. This matches section
  1: the right species' alignment is already among the candidates.
- **The knob** (`a4`; called share all / hard / FP per 20 samples): 0.5: 0.892 / 0.822 / 2; 0.35: 0.910 / 0.871 / 5;
  0.2: 0.933 / 0.960 / 8. At r226 the false positives that matter sit beside held-out congeners (v23 report). This
  world has no held-out species, so its FP column underrates the cost. Judge the knob on r226's test sets, not here.
- **Counting MAPQ < 4 pairs.** Missed strains with 2+ right-species pairs: 17 at MAPQ >= 4, 46 at any MAPQ (`a4`).
  Absent species with 2+ pairs: 10 and 141. The low-MAPQ pairs are as often a present species' reads spilling onto
  its congeners as a strain's own. As a count they would trade about 30 rescued strains for over 100 candidate false
  positives. At most, they could be a feature for a species that already has MAPQ >= 4 evidence. (An EM share over
  the candidate species is what the 2026-10-03 FN anatomy proposed; the present data do not test it.)

## 6. Memory: what more alleles or more genomes cost

**This world** ([costs.tsv](costs.tsv)). 449 species, 30 with allele genomes; the reference has 10,274,747 index
k-mers (22.9k per species). Run peak memory is 3.35-3.37 GB for every variant, dominated by the fixed key map.

| variant | alleles | edits | indexed seeds | index entries | database.protal |
|---|---|---|---|---|---|
| `reps` | 0 | 0 | 0 | +0% | 66.75 MB |
| `a4` | 10,687 | 133,270 | 0 | +0% | 66.99 MB |
| `a16` | 23,460 | 386,831 | 0 | +0% | 67.19 MB |
| `a4i01` | 10,687 (3,614 indexed) | 133,270 | 257,396 | +2.5% | 68.53 MB |
| `a4i0` | 10,687 | 133,270 | 360,575 | +3.5% | 69.14 MB |
| `a16i0` | 23,460 | 386,831 | 600,494 | +5.8% | 70.72 MB |
| `a16i0n5` | 9,180 | 144,316 | 341,718 | +3.3% | 69.02 MB |
| `oracle` (40 genomes) | 31,341 | 496,881 | 699,745 | +6.8% | 71.38 MB |

- The world's alleles look like r226's per allele. Edits per allele: 12.5 at K=4 (r226 v22: 13.2). Seeds per allele
  indexed at 0.01: 71 (r226: 76).
- **Indexed alleles roughly double the index entries of a species with 20 extra genomes** (600k seeds over 30
  species = 20k each, against 22.9k entries per species). As references of their own, the same 20 genomes would add
  20 × 22.9k = 458k entries each, 23× more.

**r226, projected** from the v22/v23 build logs: 3.12G reference k-mers; K=4 from half the genomes kept 4.22M
alleles (55.8M edits) of 7.51M offered; 1.50M indexed at 0.01 gave 114.8M seeds; the database grew by 1.1 GB and
build memory after pass 2 by 2.8 GB. Each row scales those numbers by this world's ratios. Estimates, not
measurements:

| option | alleles | index entries | memory, profiling | database |
|---|---|---|---|---|
| now (v23): K=4, half the genomes, indexed >= 0.01 | 4.2M | +115M (+3.7%) | allele table ~0.4 GB + index ~0.7 GB | +1.1 GB |
| K=4, half, not indexed (the default since `84af9a4`) | 4.2M | 0 | ~0.4 GB | ~+0.1 GB |
| K=16, all genomes, not indexed | ~10-12M | 0 | ~1.3 GB (+0.9) | ~+0.2 GB |
| K=16, all genomes, every allele indexed | ~10-12M | ~+250-300M (+8-10%) | ~+2.5-3 GB | ~+2.5-3 GB |
| up to 16 more genomes per species as references of their own | – | ~+4.8G (+155%) | ~+30 GB index, +~5 GB sequences | ~+30 GB |

How the rows were scaled:
- Alleles: K=16 keeps 98% of the offered alleles here against 45% at K=4. Doubling the genomes (the oracle's 40
  against 20) offers 1.67× as many.
- Memory of the allele table: protal's `Edit` is 4 bytes, an allele about 40.
- Index memory: packed entries of 50 bits.
- References of their own: GTDB r226 has 192,874 genomes beyond the representatives when capped at 16 per
  species (588,861 uncapped; 49,602 species have 2+ genomes, 3,805 have 17+), at about 25k entries each (215 per
  gene copy).

The build-time peak grows faster than the profiling memory: v22 needed about 24 bytes per indexed seed at build.

## 7. What to do

1. **Keep the strain alleles (scores and features).** They are the one effective lever here: +0.06 called strains,
   +0.22 for hard strains, fewer FP.
2. **Leave `--index_alleles` off.** At most a few strains in 600 are gained here, for +3.7% index entries at r226.
   The v23 build still ran with `--index-alleles 0.01`. The next build without it saves about 1.1 GB of database
   and 2.8 GB of build memory.
3. **If the final database should know more strains, use more alleles, not more references:** `--strain_alleles 16`
   and, for the final database only, all genomes (`--allele-genome-share 1`). The training database keeps the 0.5
   split for honest simulation. Measured here: +0.008 of called strains from 16 alleles instead of 4 (`a16` against
   `a4`). The share of genomes was not tested directly; 20 allele genomes instead of 5 gained +0.005. Expected cost at
   r226: about +0.9 GB of memory, unindexed. Two things need checking before switching:
   - whether a model trained on the 0.5-split database stays calibrated on a database whose species explain more
     strains (the allele features' distribution shifts);
   - what BestAllele over 16 alleles per copy costs in profiling time at r226.
4. **Do not add genomes as references of their own:** +155% index entries for at most the oracle's +0.007, and the
   MAPQ split between copies of one species would cost reads.
5. **What remains needs evidence, not references:**
   - congener-identical species (*Dorea_A*, *Priestia megaterium_D*) need reads from their discriminating sites;
   - 1-2-pair strains (the single-counted-pair group did not shrink) cannot be called without FP;
   - the 13-17 strains the model rejects with 2+ counted pairs are the ones training could reach: more real strains
     at low depth in the training samples, which the small real world had few of (23% real strains). The v23
     report's follow-up (`--insilico-strains auto`, `download_gtdb.py --per_species 4`) goes this way.

**Follow-up (same day, the user's decision).** The default became 16 alleles per species and gene (`protal --build
--strain_alleles 16`, `build_gtdb_database.py --strain-alleles 16`; 4 before). Indexing stays off by default
(`--index_alleles 1`, `--index-alleles 1`, unchanged since `84af9a4`). `scripts/ancestry_truth_test.py` still passes
its own `--strain-alleles 4`, so its documented runs stay comparable. Then the finished database's alleles came from
every genome: `build_gtdb_database.py --final-allele-genome-share 1` (new). The training database keeps
`--allele-genome-share 0.5`, and no model is trained or evaluated on the finished database. A run with no species held
out now builds a training database of its own when the shares differ; before, it trained on the finished one.

**Limits.**
- Reads from marker genes only: no off-target reads, so the samples' sizes and some sample-level features differ
  from whole-genome samples.
- One model for every variant, trained with 4 alleles from half the genomes. A model trained per variant could
  shift the `reps` and `a16` rows a little.
- No held-out species, so false positives beside novel congeners are not tested.
- Paired-end only.
- One seed: 600 strain instances per variant, so differences of a few calls are noise (see the sign tests).

**Regenerate:**
- `bash run_extract.sh` (14 min), then `bash run_all.sh` (about 30 min on 4 cores) and `bash run_oracle.sh`.
- `bash collect_costs.sh`.
- `python3 flips.py ~/strainfn <run A> <run B> ...`.

All of these run in WSL from this folder, writing to `~/strainfn`.
