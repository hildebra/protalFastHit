# Strain MSAs and trees: audit

Date: 2026-09-29. Branch `audit-fixes`, commit `7c6b2f8`. This is round 4 of the code audit that
began in [2026-09-26-code-audit.md](../2026-09-26-code-audit.md).

The audit covers how protal builds per-sample strain MSAs (`src/Profiling/Strain.h` `MSA()`,
`src/RunProtal.h` `GetMSAForTaxon`, variant calls in `src/SequenceUtils/VariantHandler.h`). It also
covers which samples and genes go into an MSA, what `scripts/qcmsa.py` filters, and how the
justfile's `strain-trees` recipe hands the result to IQ-TREE: `iqtree2 -s <species>.msa.fna -m
GTR+G -B 1000`, unpartitioned, on qcmsa's output.

## Setup

- **Build.** protal and `simulate_metagenomes` at `7c6b2f8`: a Release build in WSL Ubuntu 24.04,
  gcc 13 (`scripts/setup/setup.sh`). IQ-TREE 2.0.7 (`iqtree2`). The machine has 8 cores and is
  shared, so runs used 2–6 threads.
- **World.** `scripts/mini_db/simulate_gtdb_release.py --seed 21 --genomes_per_species 6
  --strain_divergence 0.001-0.01`, on 8 lineages:
  - 3 *Mockella* species and 2 *Fakibacter*, so congeners;
  - *Testella one* and *Dummya solo*;
  - the archaeon *Calidella fervens*.

  The database (`gtdb_to_protal_db.py`, `protal --build`) holds each species' representative and
  the shipped 2024 presence model. In this world every genome is an independent mutation of its
  species' ancestor, so strains form a star tree.
- **Known trees.** For the phylogeny work, strains were evolved from the representative's genome
  along known trees (`scripts/phylo/evolve.py`), so the reference sits at the root:
  - 12 tips per species, Yule shape;
  - K80 + Γ substitutions and small indels;
  - root-to-tip divergence 1% (*M. alpha*), 0.5% (*C. fervens*) and 0.2% (*T. one*).
- **Reads.** `simulate_metagenomes` with ART HS25, 150 bp reads, 350 ± 50 bp fragments unless
  stated. Manifests controlled every sample's strains and depths.

The work was split between four reviewers, each with its own scripts in `scripts/`:

| Part | Scripts | What was run |
|---|---|---|
| MSA construction | `msa_code/` | gtest probes against the source (`probe.cpp`); edited smoke SAMs profiled with `--profile_only` (`sam_edit.py`); cells against the truth (`truth_compare.py`, `n_reasons.py`, `overlap_depth.py`) |
| Selection, outputs, qcmsa | `qc/` | runs of 6–10 samples at 50k–250k pairs (`simrun.sh`); every qcmsa default switched off in turn (`ablate.py`); file consistency (`consistency.py`, `cellcheck.py`, `partcheck.py`); IQ-TREE as the justfile runs it (`treeeval.py`) |
| Genotype accuracy | `accuracy/` | run A: 42 samples, each genome of each species at 1, 2, 3, 5, 10, 20 and 50x, one strain per species per sample; run B: 12 two-strain mixtures (minor 2–50%, 20x or 40x); run C: 200 and 600 bp fragments; 13 reruns with filters relaxed (`run_variants.sh`); truth from affine global alignments of every strain gene to the reference gene (`truth.py`); `evaluate.py`, `compare.py`, `summarize.py` |
| Phylogeny | `phylo/` | depth series (2, 3, 5, 10, 20, 50x, uneven 1.5–100x); a two-strain sample; abundant congeners; filter and IQ-TREE variants (`analyze.py`); RF, quartet and patristic comparisons against the true trees (`treecmp.py`) |

The data (about 11 GB of reads, SAMs and runs) stays in WSL under `~/audit5` and is not in git.
Its scripts regenerate it: `setup/setup.sh` builds the world, then each part's `run_sim.sh` or
`queue*.sh` runs its simulations. qcmsa was traced with a copy of `scripts/qcmsa.py` that also
records the raw column of every output column (`make_qcmsa_colmap.py`, `make_qcmsa_trace.py`).
Its output is byte-identical to protal's.

## Result

protal's raw MSAs are accurate, and IQ-TREE recovers the true trees from them. The problems are the
qcmsa defaults and the variant calls at low depth.

- **Filtered MSA:** qcmsa's default filtering destroys strain branch lengths and drops the best
  samples. Species found in 2–3 samples get no filtered MSA at all.
- **Low depth:** the variant calls bias rows towards the reference, through the strand filter and
  through reads with no base at a position.

### Verified to hold

- **Called bases are correct:** 99.999% accurate at 10x and above, with SNP recall 0.97–0.99. On
  the smoke run, 2.68M cells held 1 false SNP. Every IUPAC cell contains the true base.
- **Rows and partitions are consistent:** all rows have the same length, and partitions are exact.
  The reference row equals the database gene, and qcmsa rewrites partitions correctly.
- **Raw MSA to IQ-TREE works:**
  - RF 0 at every depth of 3x and above;
  - tree length 0.94–1.00× the truth;
  - 99–100 bootstrap support on true splits, and no false splits.

  Partitioned trees are identical to unpartitioned GTR+G. IQ-TREE reads IUPAC codes as
  ambiguities, and including the reference row neither helps nor hurts.
- **Harmless defaults:** `--msa_min_hcov 1000` and qcmsa's `--reapply-hcov` and `--gene-min-hcov`
  never removed a row; a 1.5x row still had about 20k valid bases. Reruns give identical MSAs.

### Findings

| # | Severity | Area | Finding | Evidence |
|---|---|---|---|---|
| 1 | High | qcmsa | `--min-parsimony-samples 2` (the default) drops singleton sites, i.e. every strain's private SNPs. The reference row counts as a sample. | Terminal branches come out at 0–9% of their true length in every scenario; the raw MSA gives 87–107%. Distinct strains look identical. At 2x, 3 of 9 trees are wrong, against 1 of 9 from the raw MSA and 0 of 9 with singletons kept. |
| 2 | High | qcmsa | The multi-allelic (MRate2) filter removes clean deep samples and genes, and misses real mixtures. | Across 4 runs, 17 samples were removed; 16 were single-strain, 14 of them the deepest of their species. At 50x, 21 of 44 rows were removed. It caught 1 of 10 real two-strain mixtures. With `--snp_min_af 0` and `--snp_max_alleles 3`, errors and cross-mapped reads become about 180 IUPAC codes per clean 50x row. A pooled rate (Σ multi_allelic / Σ counts_vcov2) separates the groups: 0.45–1.88% for mixtures, at most 0.34% for clean rows. When the IQR is 0, the Tukey fence collapses to Q3. |
| 3 | High | calling | The strand filter exempts the reference allele, so a true SNP seen on one strand becomes N while the reference bases around it are kept. | At 1–5x, 13–39% of covered true SNPs become N, against 0.04–0.36% of invariant sites. This accounts for 82% of the N cells in single-strain rows. |
| 4 | High | qcmsa | `--gene-min-samples 3` is a strict `>` test, so qcmsa needs 4 samples while protal writes MSAs from 2. Species in 2–3 samples get no filtered MSA, and nothing reports it. | One species per run lost its filtered MSA, including one at 41–70x. The justfile skips the species without a message. |
| 5 | High | calling | A read rejected mid-CIGAR keeps its coverage, and any variants recorded before the rejection. | A SAM that writes mismatches as `M` (as other aligners do): 0 of 1,237 true SNPs called, 1,014 of them written as the reference base. |
| 6 | High (per species) | selection | MSA rows are admitted by the presence call, so a species the model scores near 0.5 enters a random subset of its samples. | The archaeon scores 0.42–0.54 at every depth from 2x to 50x. It entered 3 of 12 samples at 20x and 10 of 12 at 50x. The shipped model is the cause (round 3). |
| 7 | Medium | calling | Reads of relatives that map to a gene cause most false SNPs and IUPAC codes. | Cross-mapped reads are behind 97% of false SNPs and 86% of IUPAC codes in single-strain rows. Their median identity is 0.927, against 0.987 for the species' own reads. |
| 8 | Medium | calling | Overlapping mates are counted twice, in depth, observations and strands. | A single molecule passes both `--snp_min_cov 2` and the strand filter. With 200 bp fragments, 8–10% of positions with at least 2 reads come from one molecule. |
| 9 | Medium | calling | Reads with no base at a position (an N, or inside a deletion) count as reference support. | A SNP masked by N in every read is written as the reference base. Deletion-spanning reads turn a clean SNP into an IUPAC code. |
| 10 | Medium | calling | The consensus is ranked by quality sum, with the reference allele fixed at Q40 per read. | With `--snp_max_alleles 1`, the reference wins 60% of 50/50 sites, and 7 of 10 reads carrying a Q15 SNP lose to 3 reference reads. The help text says alleles are ranked by observations. |
| 11 | Medium | selection | Genes without unique k-mers are left out of the MSA. | *Dummya solo*: 59 of 118 genes, although the others have reads. |
| 12 | Medium | calling | Mismatches at read ends next to an insertion are called as SNPs. | 16% of true insertions have wrong bases or IUPAC codes next to them. |
| 13 | Medium | coverage | Reference calls need 2 reads, and a position with 1 read is written `-`, the same as a deletion. | 42% of cells are `-` at 2x. A relaxed caller at 2x (`--snp_min_cov 1 --snp_no_strand`) brings all 27 true splits to ≥95 bootstrap (from 20 of 27), but adds 116–616 false SNPs. |
| 14 | Medium | outputs | Stale filtered MSAs survive a rerun, and `strain-trees` builds trees from them. | Seen with `--msa_species` naming a species that passes nowhere, and with a species below 2 samples in the rerun. |
| 15 | Medium | qcmsa | The coverage gate counts depth ≥ 1, while the MSA writes `-` below 2 reads. The depth gate can never fire. | 187 cells passed hcov ≥ 0.3 with under 30% of the gene written. 0 of 14,129 cells were below the depth gate. |
| 16 | Medium | input | Duplicate `#SAMPLEID`s are accepted, and qcmsa then merges the two rows. IQ-TREE cuts names at their first space. | Two rows named `sample_2`, differing at 74k columns, came out with the same sequence. |
| 17 | Medium | qcmsa | `--discard-constant` counts IUPAC codes as alleles, so `+ASC` fails. | Without `+ASC`, branch lengths are inflated 38–258×. |
| 18 | Low | various | Assorted small issues, listed below. | |

The smaller issues (18):

- an insertion and a SNP at the same position compete for one call;
- `multi_allelic` in `.meta.tsv` differs from what the MSA holds;
- IUPAC ties are broken arbitrarily;
- qcmsa writes no summary on its early exits;
- the qcmsa exit status is printed as a raw value (512);
- `.meta.tsv`'s `filtered` is not the MSA's N count;
- the misc tables have no header;
- the reference row is undocumented;
- justfile: `-T AUTO` and no `-seed` in `strain-trees`; `strain-refilter` lacks `--reapply-hcov`;
  `strain-protal` uses `protal profile`, which isn't a subcommand;
- the progress bar floods log files.

### What the defaults remove

| Default | Effect on these runs |
|---|---|
| `--msa_min_hcov 1000` (rows), qcmsa `--reapply-hcov` | nothing |
| `--gene-min-hcov 0.3` | 0–1 genes |
| `--gene-min-mean-depth 1.0` | nothing; it can never fire |
| `--gene-min-samples 3` | every species found in 2–3 samples |
| MRate2 fences and cell masking | 3–6 samples and 15–64 genes per run; half the rows at 50x |
| `--min-parsimony-samples 2` | 35–70% of variable sites; terminal branches shrink to 0–9% |

### Genotype accuracy by depth (run A, bacteria except *Dummya*)

| depth | cells called, raw / qcmsa | SNP recall, raw / qcmsa | N at covered true SNPs vs invariant | false alt, ppm |
|---|---|---|---|---|
| 1 | .268 / .227 | .159 / .138 | .385 vs .0036 | 1718 |
| 2 | .577 / .508 | .390 / .346 | .322 vs .0021 | 260 |
| 3 | .783 / .692 | .594 / .530 | .239 vs .0012 | 70 |
| 5 | .945 / .836 | .821 / .732 | .128 vs .0004 | 14 |
| 10 | .994 / .880 | .969 / .864 | .019 vs 0 | 2 |
| 20 | .995 / .761 | .985 / .745 | .003 vs 0 | 1 |
| 50 | .996 / .357 | .987 / .360 | .002 vs 0 | 0 |

## Fix plan

Steps agreed on 2026-09-30, one commit each on branch `strain-fixes`:

| Step | Scope | Findings |
|---|---|---|
| L | qcmsa defaults and the tree hand-off: `--min-parsimony-samples 0`, not counting the reference row; a gene-samples floor of 2 samples; constant-site detection on A/C/G/T only; stale outputs removed; duplicate and whitespace sample IDs rejected; qcmsa summaries on early exits; justfile seed, threads and messages. | 1, 4, 14, 16, 17, parts of 18 |
| M | Variant calls for the MSA: strand test against the strands of all reads at the site (reference included); each fragment counted once; depth from reads with a base; consensus by observations; rejected reads leave nothing; read-end mismatches next to indels ignored; insertions kept apart from SNPs. | 3, 5, 8, 9, 10, 12 |
| N | Noise and relatives: an Illumina `--snp_min_af` for IUPAC codes (this leaves profiles unchanged, since the model's features ignore variant validity); an identity margin for reads used in variant calls; MRate2 on the pooled rate, without iteration, guarding IQR = 0. | 2, 7 |
| O | Coverage and genes: positions with one read are called; genes chosen by their reads; qcmsa's gates on what the MSA holds. | 11, 13, 15 |
| P | Which samples enter an MSA: their own `--msa_knob`, defaulting to `--knob`; species whose own reads give strong evidence but whose score stays below the knob are reported. | 6 |
