# The depth identity margin scaled per gene, and congeners

- **Date**: 2026-10-01.
- **Question**: the [stress test of 2026-09-30](../2026-09-30-depth-margin-stress/README.md) found that, where genes
  differ in conservation, the fixed margin (98th percentile − 0.08) undercounts strains whose reads on fast genes
  fall further below the top. (1) Scale the default 0.08 per gene by the gene's conservation, and compare it with
  the gene-scaled median; (2) test both on a world with many congeners, many of them missing from the database;
  (3) two rules for relatives' reads raised on the way: a narrower margin for species whose congeners were
  detected, and a check of a species' depth on conserved against fast genes.
- **Code**: `audit-fixes` at `25354fd` (this change: the factors and the scaled margin); the other rules are an
  experimental build of it (`scripts/rule_patch.py` patches a copy of the source; the rule from
  `$PROTAL_DEPTH_RULE`) that profiles the same SAMs again (`--profile_only`).
- **Machine**: WSL2, Intel Core Ultra 7 258V, 6 threads, 23 GB; shared during the runs with another session's
  benchmarks (the runs here were niced; timings in this report are of builds only).
- **Data**: three simulated worlds (`simulate_gtdb_release.py`, 5 genomes per species, 200 kb genomes; samples
  by `simulate_metagenomes --from_manifest`, ART HS20 2x100 and HS25 2x150, 6 samples each; for each world two
  databases built by this protal, every species (`db7_full`) and without the held-out species (`db7_missing`)):

  | World | Species | Genes | Strains from the reference | Congeners | Held out | Samples |
  |---|---|---|---|---|---|---|
  | 1 (`~/stress`, 2026-09-30) | 120 GTDB-like | one rate | 0.4–6% | 3–12% apart | 30, each with a congener left | 60 species each, a quarter of each kind |
  | 2 (`~/stress2`, 2026-09-30) | the same | `--gene_rates categories` | 0.4–6% | 3–12% at a gene of factor 1 | the same | the same |
  | congener (`~/stress3`, new) | 160: 12 genera of 10, 40 in small genera | `--gene_rates categories` | 0.4–4% | 2–10% at a gene of factor 1 | 4 of every large genus (48) | 4 large genera whole and 20 other species each; held-out species 5 times as deep in half of them |

  Kinds of species in every sample, in turn: the representative alone, one other strain alone, the
  representative with the most distant strain as a 10–30% minor, and the reverse. The congener world:
  `scripts/congener_lineages.py`, `congener_world.sh`, `congener_design.py` (`results/design_congener.tsv`,
  `results/heldout_congener.txt`).
- **Run**: `scripts/run_rules.sh WORLD` (databases, simulation, alignment once per database with `--no_profile`,
  then each rule with `--knob 0 --no_strains`); `scripts/score.py WORLD` (`results/scores_*.md`);
  `scripts/factor_check.py`, `scripts/split_diag.py`. The first two worlds' SAMs of 2026-09-30 were profiled
  again (the reference is the same); the rules that use no gene factors (no filter, 0.04, gene median − 0.08)
  are their profiles of then.
- **Scores**: as in the stress test: with `--knob 0` every taxon with reads is reported, so the scores are of
  abundance, not detection. Bray-Curtis between the true relative abundances of the species present (less the
  held-out ones against `db7_missing`) and the reported ones renormalised over them; median log2(reported /
  true) by kind of species and by its congeners in the sample.

## What was implemented

- **Factors at `--build`** (`src/SequenceUtils/GeneConservation.h`, `Build.h`): for each species and gene, up
  to 16 other genomes' copies in `full_reference.fna` against the representative's gene (Mash distance of the
  12-mer sets, copies more than 0.25 away left out, one identical copy taken as the representative's own); divided
  by the distance of the species' median gene; the factor of a gene is the median over the species with 10 genes
  or more and a median gene 0.2% or more from the representative, scaled so that the median gene has 1, shrunk
  towards 1 by 3 species' weight, within 0.25–4. Written as `gene_conservation.tsv` (gene id, factor, species)
  into the database, a part of `database.protal`; the log reports it (`Gene conservation: factors 0.26-3.3 for
  168 genes, from 120 species`). Without other genomes' copies there is no table.
- **Scaled margin** (`Profiler.h`, `Taxon::OwnIdentityThreshold(geneid)`): a read on gene g counts towards the
  depth if its identity is at least the 98th percentile − (0.03 + (M − 0.03) × r_g), M =
  `--depth_identity_margin` (0.08): 0.03 for read errors and the spread of read identities, the rest scaled.
  `--gene_conservation FILE` uses other factors, `--gene_conservation none` the same margin on every gene, as
  databases without the table do.
- **Tests**: unit tests of the margin, the table's reading and its errors, the k-mer distance, the estimate on
  simulated genes of three rates, and the depth of a species with a fast and a conserved gene; e2e tests of a
  build's factors, of the queries with and without them, and of a malformed table.

The factors recover the simulated rates (`results/factors_*.txt`, true rates scaled to a median of 1):

| World | Genes | Factors | Correlation with the true rates | ribosomal / translation / tRNA / other, estimated (true) |
|---|---|---|---|---|
| 1 (one rate) | 168 | 0.87–1.1 | (true rates all 1) | — |
| 2 | 168 | 0.26–3.3 | 0.985 | 0.44 (0.38) / 0.85 (0.80) / 1.31 (1.24) / 1.30 (1.31) |
| congener | 120 | 0.25–3.2 | 0.999 | 0.50 (0.49) / 0.88 (0.88) / 1.22 (1.24) / 1.42 (1.45) |

The step took 1.3 s of a 5.6 s build of the congener world (94,063 copies, all compared: 5 genomes per
species). At GTDB r226 size it is one more pass over `full_reference.fna` (~51 GB, extrapolated from the
tuning world) with at most 16 copies compared per species and gene: an estimated 1–4 minutes of a build of
hours.

## Results

Bray-Curtis, mean of 6 samples per cell; the sum is that of the four cells
(`results/scores_congener.md`, `scores_world2.md`, `scores_world1.md`):

| Rule | congener: full 2x100 / 2x150 | congener: missing 2x100 / 2x150 | congener sum | world 2 sum | world 1 sum | all |
|---|---|---|---|---|---|---|
| no filter | 0.058 / 0.038 | 0.150 / 0.163 | 0.409 | 0.165 | 0.174 | 0.748 |
| 98th percentile − 0.04 (0.7 until 1b6161d) | 0.149 / 0.088 | 0.136 / 0.086 | 0.459 | 0.621 | 0.536 | 1.616 |
| 98th percentile − 0.08, the same on every gene | 0.061 / 0.040 | **0.133 / 0.134** | 0.368 | 0.182 | **0.162** | **0.712** |
| **scaled: 0.03 + 0.05 r (the new default)** | 0.059 / 0.038 | 0.140 / 0.153 | 0.390 | 0.169 | 0.165 | 0.724 |
| scaled: 0.08 r | 0.060 / 0.038 | 0.139 / 0.152 | 0.389 | 0.178 | 0.165 | 0.732 |
| scaled: 0.02 + 0.06 r; 0.04 + 0.04 r | 0.060 / 0.038; 0.059 / 0.038 | 0.140 / 0.153; 0.141 / 0.154 | 0.391; 0.391 | | | |
| scaled: 0.03 + 0.07 r (margin 0.10) | 0.058 / 0.038 | 0.148 / 0.161 | 0.405 | | | |
| gene median − 0.08 | 0.058 / 0.038 | 0.149 / 0.162 | 0.407 | 0.161 | 0.169 | 0.737 |
| gene-scaled median (M 0.03, K 0.03) | 0.058 / 0.038 | 0.147 / 0.160 | 0.404 | **0.159** | 0.164 | 0.727 |
| gene-scaled median (M 0.03, K 0.02) | 0.059 / 0.038 | 0.145 / 0.158 | 0.400 | | | |
| detected congener → margin 0.05 | 0.083 / 0.050 | **0.103 / 0.086** | **0.322** | 0.283 | 0.273 | 0.878 |
| detected, deeper congener → margin 0.05 | 0.075 / 0.047 | 0.119 / 0.135 | 0.376 | 0.238 | 0.237 | 0.851 |
| conserved genes alone if < 0.8 × the fast ones' depth | 0.109 / 0.081 | 0.111 / 0.112 | 0.412 | 0.211 | 0.165 | 0.788 |
| the same at < 0.7 | 0.102 / 0.072 | 0.104 / 0.111 | 0.389 | 0.189 | 0.165 | 0.743 |

In the congener world, by the species' congeners in the sample (median log2(reported / true) / mean |log2|):

| Rule | full: closest congener < 4% (n=239) | missing: held-out congener present (n=288) | missing: outnumbered by one (n=194) | full: strain alone 2–4% (median) |
|---|---|---|---|---|
| no filter | −0.030 / 0.225 | +0.116 / 0.553 | +0.408 / 0.729 | −0.121 |
| 98th percentile − 0.08 | −0.030 / 0.234 | +0.102 / 0.480 | +0.266 / 0.628 | −0.123 |
| scaled (new default) | −0.029 / 0.228 | +0.100 / 0.519 | +0.314 / 0.682 | −0.124 |
| gene-scaled median | −0.035 / 0.223 | +0.084 / 0.546 | +0.355 / 0.724 | −0.121 |
| detected congener → 0.05 | −0.060 / 0.293 | +0.056 / 0.358 | +0.127 / 0.440 | −0.242 |
| conserved genes alone (< 0.7) | −0.087 / 0.437 | +0.097 / 0.411 | +0.268 / 0.514 | −0.076 |

- **The scaled margin helps where strains are far, and loses where relatives are missing.** In world 2, whose
  strains reach 6% from the reference, it scores 0.169 against 0.182 for the same margin on every gene (strains
  4% or more away: −0.228 against −0.242 median log2). In the congener world, whose strains are within 4% but
  whose missing congeners are 2–10% away and often dominant, it scores 0.140 / 0.153 against 0.133 / 0.134 with
  species missing: a wider margin on fast genes admits the missing relative's reads there too, since its reads
  sit below the top in proportion to the gene's rate just as a strain's do. With every species in the database,
  it is as good as no filter (0.059 / 0.038).
- **Where genes evolve alike** (world 1), the build's factors are 0.89–1.1, noise around the true 1, and the
  scaled margin costs 0.003 (0.165 against 0.162): what the estimate's noise costs when there is nothing to
  scale.
- **How the margin is split barely matters.** The fixed part 0 to 0.04 of 0.08 scores 0.389–0.391 in the
  congener world; 0.08 r is worse than 0.03 + 0.05 r in world 2 (0.178 against 0.169).
- **The gene-scaled median** is the best rule in world 2 (0.159) and close to no filter in the congener world
  with species missing (0.147 / 0.160): a relative that outnumbers a species pulls the median over genes with
  it, as the first world found.
- **A narrower margin for species with a detected congener** (the model's score 0.5 or more, under the default
  rule; a two-pass rule) is the best rule in the congener world (0.322), where almost every species of a large
  genus has one: the missing congeners' reads go (+0.127 against +0.266 for the species they outnumber). It costs
  strains: 2–4% from the reference −0.242 against −0.123 in the congener world, and in worlds 1 and 2, whose
  strains reach 6%, it is the worst rule but 0.04 (0.273, 0.283). Narrowing only for a deeper detected congener
  costs less and gains less (0.376, 0.238, 0.237). Over the three worlds it is fitted to the one it wins.
- **The conserved and fast genes' depths: the premise was wrong.** A relative the database lacks does not make a
  species' conserved genes deeper, but its fast ones (`results/split_diag_*.txt`): the ratio of the conserved
  genes' median depth to the fast ones' (all reads) is 0.71 for species outnumbered by a held-out congener,
  0.84 beside one that does not outnumber them and 0.93 without one (congener world, `db7_missing`; world 2:
  0.89, 0.92, 0.97). Presumably (not measured) the relative's reads on a conserved gene match several congeners
  about equally and are dropped as ambiguous (MAPQ), or land on genes without unique k-mers, while on a fast gene
  they align uniquely to the nearest congener in the database. The stress test's statement that "a relative's reads now align mostly on the
  conserved genes" was reasoned, not measured, and is wrong. Turned round (the conserved genes alone when they
  are less than 0.7–0.8 times as deep as the fast ones), the check halves part of the error with species missing
  but fires as often without a missing relative: congeners in the database lower the conserved genes too (0.86
  in `db7_full`), and its sums are no better than the fixed margin's. *Note (later on 2026-10-01): these ratios
  are of present species after the profiler's filters; they do not say where a relative's reads align. Traced read
  by read on the v0.7.1 benchmark world, a missing relative's reads align best on the conserved genes, most of them
  there fitting several congeners equally (MAPQ below 4)
  ([conservation pattern report](../2026-10-01-conservation-pattern/README.md)).*

## Recommendation

- **Over the three worlds, the same 0.08 on every gene is the best rule** (0.712), then the scaled margin
  (0.724) and the gene-scaled median (0.727): 1–2% apart, less than the worlds differ. Scaling the margin pays
  only where strains reach 6% from the reference and relatives are rarely missing (world 2); it loses where
  many congeners are missing (the congener world) and, by its noise, where genes evolve alike (world 1). The
  default of `25354fd` (scaled) should go back to the same margin on every gene, with the scaling kept as an
  option, unless real genomes say otherwise; that is for the team to decide.
- **Keep the factors in the database.** They cost minutes of a GTDB build, recover simulated rates at 0.985–0.999,
  and are what the scaled margin and the gene-scaled median need, as options and for a check on real genomes.
- **No congener-aware rule as a default.** The narrower margin beside detected congeners wins only the world of
  large genera with missing members and costs distant strains everywhere else; the conserved-against-fast check
  rests on a premise the data reverse, and turned round it fires as often without a missing relative.
- **What would decide it**: real genomes, a GTDB genus of many species with some held out, and genes whose rates
  differ between lineages, which none of these worlds has. The errors that remain largest are not the margin's:
  species outnumbered by a missing congener (+0.27 median log2 with the fixed margin, +0.13 with the narrower
  margin beside detected congeners), and strains 4% or more from the reference (world 2: −0.17 to −0.27 under
  every rule), whose reads are lost at alignment.

## Follow-up: the default reverted

As decided the same day, the margin is the same on every gene again: `--gene_conservation` defaults to
`none` (`25354fd` had made the scaled margin the default), and `--gene_conservation db` scales it by the
database's factors. `--build` still estimates the factors and stores them in every database built with
other genomes' copies, so the scaled margin, or any other use of the genes' rates, needs no rebuild. A
query's log says whether the database has them, `protal --unpack_db` writes the table out, and
`build_gtdb_database.py` records the build's summary of them in `build_metadata.tsv` (`gene_conservation`).
The results above stand: the default is the rule "98th percentile − 0.08, the same on every gene" again, and
the one marked "the new default" is the option.
