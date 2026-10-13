# A local test of the mock benchmark: ten real samples on a reduced r226 database

2026-10-13. protal and simulate_metagenomes at `4dfc549` (v0.8.3), built in a WSL clone. The mock benchmark's scripts
are from `scripts/mock_benchmark/` of the
[survey](../2026-10-13-public-mock-communities/README.md) (`1dcc325`, `f072328`), plus the change below.

**Machine and limits.** WSL on the user's laptop: 6 vCPUs, 23 GB of memory. Everything ran on 4 cores
(`taskset -c 0-3`). The full r226 database (~35 GB in memory) does not fit, so this test builds a reduced
one, on the user's instruction: from the local GTDB r226 folder, holding the species representatives of the mock
communities.

**What it tests.**
- **The pipeline end to end on real reads:** download, a database, `score_mocks.py run` and `score`.
- **What protal makes of these reads:** the scores, the errors behind them, and which errors are protal's.
- **Not protal at r226.** The reduced database lacks most of GTDB, so false positives from absent relatives are
  undercounted, and its models are trained on a small world.

## Setup

**Database: 761 GTDB r226 species** ([select_mock_taxa.py](select_mock_taxa.py), [make_inputs.sh](make_inputs.sh)).
- **Species:**
  - **Truth species:** all 269 species of the 29 truth tables.
  - **Congeners:** up to 3 more species of each of their genera, as competitors (176 genera).
  - **Strains:** up to 5 strain genomes per species (CheckM2 ≥ 90/≤ 5, isolates first).
- **Inputs:**
  - **Marker genes:** the species' marker genes come from the local r226 archives (`extract_subset.py`); 2,634 genomes
    in all (bac120 2,475, ar53 159).
  - **Genomes for the training simulations:** from NCBI (`fetch_genomes.py`): the representatives of 761 species and
    277 strains of 150 species. 1,028 of 1,038 arrived; 10 UHGG-derived MAG accessions are no longer served.
- **Build** ([run_build.sh](run_build.sh)): `build_gtdb_database.py` scaled to 4 cores, as the small real world's
  ([2026-10-10](../2026-10-10-real-small-world/README.md)), with Nanopore reads added.
  - **Training samples:** 200 pe, 48 HiFi and 48 Nanopore.
  - **Test samples:** 30 pe, 16 HiFi and 16 Nanopore.
  - **Hold-outs:** 30% of the species and 3 genera.
  - **Gates:** the default sparse-first gates (pe 4, pb 16, ont 8), with `--model-gates none`.
  - **Result:** 14 min 22 s, peak 4.1 GB, `database.protal` 151 MB ([build_steps.txt](results/build_steps.txt)).
  - **Its own test sets** ([build_model_summary.txt](results/build_model_summary.txt)):

    | Read type | F1 | False positives per sample |
    |---|---|---|
    | pe | 0.963 | 0.47 |
    | HiFi | 0.968 | 0.44 |
    | ONT | 0.962 | 0.25 |

**Samples.** Ten samples were drawn at random by [pick_samples.py](pick_samples.py) (seed 20261013): ten of the 12 sets
that have a truth table, and one sample of each.
- **Download** ([run_test.sh](run_test.sh)): the first 1M read pairs of each short-read sample and the first 1 Gb of
  each long-read sample, 3.8 GB in all, into `~/Mocks`.
- **Two more paired-end D6331 samples** of other labs (see [Dohlman sample](#the-dohlman-d6331-sample-is-an-outlier)),
  added by [run_more_pe.sh](run_more_pe.sh) with [samples_extra.tsv](samples_extra.tsv):
  - HKU's SRR23926922;
  - Auckland's SRR33710519.
- **Runs:** protal ran on 4 threads into the run folder `~/Mocks/0.8.3_4dfc549`, named from its `--version`. Each
  sample took 1-32 s at a peak of 3.3-3.4 GB.
- **Ungated rerun** ([run_ungated.sh](run_ungated.sh)): the same ten samples again with `--sparse_first 0`, into
  `0.8.3_4dfc549_ungated` (`--run-name`, `--protal-args`).
- **Order:** [run_all.sh](run_all.sh) runs these steps in order.
- **Tables:** in [results/](results/): `samples_by_run.tsv`, `summary_by_run.tsv`, and `species_0.8.3_4dfc549.tsv`
  (truth and profile shares per species).

## Results

Run `0.8.3_4dfc549`; the ungated run is the same within noise (see [Not the gate](#not-the-gate)). Shares compare cells,
over the species the truth places.

| Sample | Reads | Truth (kind) | Called | TP | FP | FN | F1 | L1 | Bray-Curtis | `?` (expected) |
|---|---|---|---|---|---|---|---|---|---|---|
| D6331 pe, Dohlman (DNA only) | 1M pairs | abundance | 11 | 11 | 0 | 4 | 0.85 | 0.92 | 0.46 | 0.45 (0.006) |
| D6331 pe, HKU (extra) | 1M pairs | abundance | 12 | 12 | 0 | 3 | 0.89 | **0.19** | **0.09** | 0.13 (0.006) |
| D6331 pe, Auckland (extra) | 1M pairs | abundance | 12 | 12 | 0 | 3 | 0.89 | 0.41 | 0.20 | 0.20 (0.006) |
| D6331 in human DNA 1:10^8 | 1M pairs | abundance | 0 | 0 | 0 | 15 | 0 | 1 | 1 | 1 |
| Treichel Mock-even-70 + mouse DNA | 1M pairs | abundance | 63 | 62 | 1 | 1 | 0.98 | 0.46 | 0.23 | 0.76 (0) |
| D6331 HiFi, ultra-low input (WGA) | 1 Gb | abundance | 15 | 13 | 2 | 2 | 0.87 | 0.77 | 0.38 | 0.14 |
| ATCC MSA-1003 HiFi | 1 Gb | abundance | 18 | 17 | 1 | 3 | 0.89 | 0.46 | 0.23 | 0.27 (0) |
| D6331 ONT R10.4.1 sup v5.2 | 1 Gb | abundance | 14 | 14 | 0 | 1 | 0.97 | 0.47 | 0.23 | 0.15 |
| 3 *E. coli* strains, ONT R10 | 0.67 Gb | closed | 1 | 1 | 0 | – | – | – | – | 0.11 |
| NIBSC multi-lab negative control | 40 pairs | negative | 0 | – | 0 | – | – | – | – | 1 |
| hCom2-like in vitro, passage 2 | 1M pairs | closed (117) | 52 | 52 | 0 | (65) | – | – | – | 0.11 |
| hCom1a mouse caecum | 1M pairs | dropout | 54 | – | – | – | – | – | – | 0.19 |

**How to read the table.**
- **Not real misses:**
  - **D6331's three rarest species:** *S. enterica*, *E. faecalis* and *C. perfringens* (7e-5 to 9e-7 of cells) are
    undetectable at these depths.
  - **The human-DNA sample at 1:10^8:** it holds no microbial reads at 1M pairs, so its "15 FN" are expected.
  - **The strain mix's "FN":** it is the closed set's other species, *K. pneumoniae*, which this mix does not contain.
- **The hCom2 closed set:** a member missing at passage 2 is not an error (the truth is the membership only). What
  matters: 52 members called, no species outside the set.
- **hCom1a:** both strains the dropout arms remove, *C. scindens* and *E. hylemonae*, are called in the full community.
- **False positives:** only three species, all minor or congeners:
  - *Bifidobacterium catenulatum* at 0.02% in Treichel's mock;
  - `Escherichia sp004211955` in both HiFi samples (below);
  - a 0.01% *Veillonella* congener in the D6331 HiFi sample.
- **Controls:** the blank and the human-DNA sample have no calls.

### *E. coli* missed in HiFi reads, a congener called instead

In both HiFi samples *E. coli* is not called, and `s__Escherichia sp004211955` is called in its place:

| Sample | *E. coli* depth | Probability | Reads' identity | *Escherichia sp004211955* |
|---|---|---|---|---|
| ATCC MSA-1003 (*E. coli* K-12 MG1655, 10% of cells) | 46.5× | 0.122229 | 98.2% | 2.5×, P 0.55 |
| D6331 HiFi (five *E. coli* strains) | 25.1× | 0.122229 | 98.2% | 1.2×, P 0.63 |

Both samples get the same probability to six digits, so the gradient-boosted model's decision for *E. coli* ends in
the same leaves in both: their features fall into the same bins.

- **ATCC:** *E. coli* is the sample's most abundant species, with reads on 81 of 103 marker genes and a mean MAPQ of 15.
- **The cause:** the representative of r226's *E. coli* is GCF_003697165.2, not K-12. The reduced database holds only
  5 more *E. coli* strains, none of them K-12. So K-12's reads are 98.2% identical to the genes they align to, and
  the HiFi model takes a species at that distance for a relative.
- **The other read types call it:** the ONT model calls *E. coli* in D6331 (P 0.59, identity 97.3%), and the pe model
  calls it in the HKU and Auckland samples (P 0.80 and 0.82, 97.7-97.9%).

This is the error class of divergent strains seen in the r226 evaluations. Whether it happens at r226, whose
*E. coli* has thousands of genomes and up to 16 strain alleles per species, is the first thing to check on the
cluster with the full database. The case is clean: any HiFi sample of ATCC MSA-1003 or Zymo D6331.

### The Dohlman D6331 sample is an outlier

The Dohlman sample alone has four species far below their share, with reads of lower identity:

| Species | Dohlman: depth relative to *B. fragilis* | HKU | Auckland | Truth (cells) |
|---|---|---|---|---|
| *V. rogosae* | 0.05 (identity 97.0%) | 1.40 (98.1%) | 1.60 (98.3%) | 2.4 |
| *R. hominis* | 0.05 | 0.86 | 0.28 | 1.5 |
| *P. corporis* | 0.03 | 0.72 | 0.76 | 0.75 |
| *E. coli* | 0.14, not called (P 0.016) | 0.64 | 0.51 | 1.05 |

- **Same pattern in the long reads:** the same standard on HiFi and on ONT gives *V. rogosae* about twice
  *B. fragilis*, as the truth does.
- **The cause is the sample, not protal's paired-end path:** the same database and binary profile the HKU sample
  within L1 0.19 of Zymo's manual (*E. coli* 8.5% for 8.7%, *V. rogosae* 18.4% for 19.9%).
- **Likely explanation:** another lot or other strains (lower identity for exactly these four species), which the
  study's data alone cannot settle.
- **For the benchmark:** the `d6331_dna` and `d6331_host` sets should not be scored for abundance against the manual
  until this is understood. Their detection limits (the dilution series) are still usable.

### Not the gate

The ungated run (`--sparse_first 0`) gives the same calls on 9 of the 10 samples, and L1 within 0.0015 on all of them.
The one difference is hCom2, with 50 calls instead of 52. The default gate of `4dfc549` (pe 4, stride pick) does not
lose these samples' reads, divergent strains included.

### Abundances

- **Best case:** the HKU pe sample is close to the manual (Bray-Curtis 0.09). HKU and Auckland, the same standard from
  two labs, differ by Bray-Curtis 0.16 with identical calls.
- **Elsewhere, L1 0.46-0.77.** Several things add up:
  - extraction or lot effects in the cell standards;
  - the missed *E. coli* (10% of the ATCC cells);
  - the amplified library of the D6331 HiFi sample (WGA, L1 0.77);
  - *Cereibacter sphaeroides* at 28% instead of 11% in the ATCC DNA mock. It has two chromosomes and duplicated
    marker genes, as noted in the survey.
- **Unknown share (`?`):** above its expected value in every D6331 sample (0.13-0.20 for the good pe samples, against
  0.006 of yeast cells). On the build's simulated pe test samples it is within 0.02 of the truth (median -0.019), so
  this is a real-read effect. Not resolved here; plausible causes:
  - reads of real strains that align less often than simulated ones, so that the depths are underestimated and their
    bases left over;
  - the reads protal saw of species it did not call.

  Treichel's 0.76 includes the germ-free mouse DNA added to that library.

## Fixed in the scripts along the way

`score_mocks.py run` looked for each sample's profile in `profiles/` only. But a map without `SAM` and `PROFILE` columns
makes protal write `<PREFIX>.sam.zst` and `<PREFIX>.profile` into `#OUTPUT_DIR` itself. That is as `--map_help` says,
while it also creates `alignments/` and `profiles/` and leaves them empty, and docs/running.md's table reads as if the
profiles went to `profiles/`. The fix:
- **Map columns:** the script now writes the map with `SAM` and `PROFILE` columns, so the files go to `alignments/` and
  `profiles/`.
- **Finding profiles:** it finds a run's profiles anywhere under its folder.
- **Tests:** they cover both layouts.

## Next

1. **ATCC MSA-1003 and Zymo D6331 HiFi at r226 on the cluster:** whether *E. coli* is called there (see above).
2. **The Dohlman set:** score `d6331_dna` and `d6331_host` for detection only, and use the HKU and Auckland D6331 runs
   for abundance. Both could join `samples.tsv` as a set of their own.
3. **The full starting sets** on the cluster with the r226 database, as planned in the survey.

**Reproducing it.** `bash run_all.sh 4dfc5496c976032b429b081a04406411569b961b` in WSL, then `run_ungated.sh` and
`run_more_pe.sh` (scratch copies of `samples.tsv`).
- **Kept until the user decides:**
  - `~/mocktest`: the clone (0.26 GB), the inputs (1.9 GB) and the build (2.0 GB with its genome store);
  - `~/Mocks`: the reads and both runs (4.0 GB).
