# Mini database generation: which genes, and how they evolve

Date: 2026-10-05. Code at `00bebbf` (audit-fixes): `scripts/mini_db/` (`build_mini_db.sh`,
`simulate_gtdb_release.py`, `markers_r226.tsv`, `gtdb_to_protal_db.py`), `scripts/rank_genes.py`. Real-data
comparison: `local/v10/model_logs/gene_congeners.tsv` of the r226 v10 build (each real gene's conservation factor),
gene ids mapped through a mini database's `gene2geneid.tsv` (WSL `~/protal_db_v1`, same marker table, same
converter order). Script: [`scripts/rates_vs_r226.py`](scripts/rates_vs_r226.py), output
[`rates_vs_r226.tsv`](rates_vs_r226.tsv).

## How a mini database is made

`build_mini_db.sh` runs four steps:

1. `simulate_gtdb_release.py` writes a synthetic release in GTDB's layout (taxonomy, metadata, per-marker FASTAs of
   representatives and of all genomes, genome FASTAs) plus truth files (`simulation/`). Default: 3 species, 3
   genomes each, 150 kb of background DNA in 3 contigs.
2. `gtdb_to_protal_db.py` converts it as it converts a real release.
3. `gene_neighbours.py` places the genes in every genome.
4. `protal --build`.

The benchmark and training worlds use the same generator with `gtdb_like_lineages.py` (900 species), ranges of
divergence and, in some, `--operons` and `--gene_rates categories`:

| world | options |
|---|---|
| `build_mini_db.sh` default (e2e, CI-style checks) | 3 species, strain 0.005, species 0.035, all genes alike, shuffled order |
| model-training tuning world (2026-09-29) | 900 species, strain 0.002-0.02, species 0.015-0.06, all genes alike |
| operon world (2026-10-02 neighbours) | as tuning + `--operons` |
| 0.7.x benchmark world (2026-10-01 v071, 0.7.5 bench) | as tuning + `--operons --gene_rates categories` |
| gene-scaled margin / depth-margin stress worlds | `--gene_rates categories` |

## Which genes

**The marker universe is fixed** by `markers_r226.tsv`: GTDB r226's `bac120` (120) and `ar53` (53) marker
sets with their HMM lengths, 173 rows, 168 distinct markers (5 are in both sets: `PF00410.20` S8, `PF00466.21`
L10, `TIGR01171` rplB, `TIGR00064` ftsY, `TIGR00967` secY). Every species carries the full set of its domain:
bacteria the 120, archaea the 53. There is no selection, apart from:

- `--marker_loss` (default 0.02): each genome loses each marker independently with probability 2%, the same for
  every gene and clade.
- The converter assigns gene ids by first appearance: the bac120 marker files sorted by id, then the ar53-only
  ones. The 5 shared markers get one id, so a bacterium's and an archaeon's S8 are the same protal gene, as in a
  real build. The mini and real r226 databases therefore have the same 168 gene ids in the same order, which is
  what made the comparison below possible.

**Sequences.** Each marker starts as one random coding sequence (ATG, random sense codons, stop) of 1.0-1.15x
its HMM length. It evolves down the lineage (per-site rates: 0.15 into a domain, 0.10 phylum, 0.04 class, 0.03
order/family/genus, `--species_divergence` into a species, `--strain_divergence` per genome). Substitutions never
create a stop codon and 2% of events are codon indels. There is no codon-position bias: all three positions
change alike, unlike real genes and unlike the in-silico strains of `insilico_strains.py`.

**Gene speeds.** By default every gene evolves at the same rate, so all conservation factors are ~1 and the
conservation features carry no gene-to-gene signal. `--gene_rates categories` multiplies each gene's rate by a
category picked from its name (ribosomal 0.4, translation/transcription 0.8, tRNA synthetases and modification
1.1, other 1.4), times gamma noise (CV 0.35), normalised to a mean of 1.

## Findings

### 1. The categories rank genes in the right direction but get the speeds wrong

Against the real r226 factors (168 genes): Spearman +0.52 with the within-species factor, +0.50 with the
between-congener factor.

| category | genes | mini rate | real r226 factor, median (quartiles) |
|---|---|---|---|
| ribosomal | 36 | 0.4 | 0.48 (0.39-0.72) |
| translation | 31 | 0.8 | 0.96 (0.77-1.09) |
| tRNA | 21 | 1.1 | 1.34 (1.27-1.40) |
| other | 80 | 1.4 | 1.04 (0.93-1.26) |

The real order is ribosomal < translation ≈ other < tRNA. The mini world makes "other" the fastest; in reality
the tRNA synthetases and modification enzymes are. The real factors span 0.28-1.57; the mini rates, the category times gamma noise, span roughly 0.2-2.

### 2. Twelve archaeal ribosomal proteins are not recognised as ribosomal

`marker_category()` matches names by prefix (`ribosomal_`, `rps`, `rpl`, `us11`, ...). GTDB's ar53 TIGRFAM models
use the newer universal names, which none of the prefixes catch:

- **Ribosomal proteins categorised "other"** (rate 1.4 instead of 0.4, 3.5x too fast): `arch_S11P`,
  `uL14_arch`, `uS9_arch`, `uS13_arch`, `uS10_euk_arch`, `uS12_E_A`, `uS2_euk_arch`, `uL16_euk_arch`,
  `uS5_euk_arch`, `uS4_arch`, `uS7_euk_arch`, `uS3_euk_arch`. Their real factors are 0.68-0.99, below their
  set's median gene.
- **Archaeal translation and transcription factors categorised "other"**: `EF-1_alpha`, `eif2g_arch`,
  `EIF_2_alpha`, `aIF-2`, `eIF_5A`, `eIF-6`, `aRF1/eRF1`, `pelota`, `RNA_pol_A_bac`, `RNA_pol_rpoA1/A2`,
  `KOW_elon_Spt5`.
- **Bacterial misses:**
  - `secE_bact` is "translation" (0.8) but is one of the most conserved genes (0.40).
  - The ATP synthase genes (`atpD`, `ATPsyn_F1gamma`, 0.76-0.91) are "other".

So in the worlds built with `--gene_rates categories` (the 0.7.x benchmark world among them), archaea's
conserved-gene pattern is partly scrambled: a quarter of their markers evolve as fast genes though they are among
their slowest. Any archaeal result on those worlds that rests on the conservation features
(`conserved_fast_*`, `excess_scaled_*`, the per-gene margin) is measured on a pattern real archaea do not have.
Bacteria are less affected: their ribosomal and tRNA names are caught.

### 3. A few genes are much shorter than the real ones

Gene length comes from the HMM model length. The 152 TIGRFAM models are full-length (mean 338 aa). Several of
the 16 Pfam models are domains of a longer protein:

| model | mini length | real protein |
|---|---|---|
| `PF03726.15` PNPase | 83 aa | polynucleotide phosphorylase, ~700 aa |
| `PF01000.27` RNA_pol_A_bac | 112 aa | RpoA, ~330 aa |
| `PF02576.18` DUF150 | 73 aa | RimP, ~150 aa |
| `PF07541.13` EIF_2_alpha | 114 aa | aIF2α, ~270 aa |

This assumes GTDB's marker FASTAs hold the whole called gene, as they appear to; check on the HPC with
`awk '{n[$2]++; s[$2]+=$4-$3} END{for(g in n) print g, s[g]/n[g]}' protal_db/reference.map` against
`gene2geneid.tsv`. Short genes take fewer reads and have fewer k-mers, so in the mini world these genes weigh
less in profiling and score lower in `rank_genes.py` than they would at GTDB. That is four genes of 168: a
small effect on profiles, but it matters for a reduced database's gene choice, which length enters through
prevalence and unique k-mers.

### 4. What the mini worlds can and cannot say about gene selection

`rank_genes.py` (the reduced databases, `--n-genes`) scores each gene by prevalence x unique k-mer share, with
a third of the picks reserved per domain. In a mini world:

- **Prevalence** is 98% for every gene (uniform marker loss): it never discriminates. In GTDB it does (bac120
  genes are found in different shares of species).
- **Unique k-mer share**, without `--gene_rates`, varies only with length and noise. With categories it follows
  the wrong order of finding 1, so the "most distinctive" genes would be the "other" genes rather than the tRNA
  synthetases.

So the `GeneSubsetTest`/`GtdbBuildTest` reduced-database runs test the machinery (ids, files, per-domain
reservation, archaea covered). They do not test which genes a real ranking would pick. Only a ranking from a real
build (`build_gtdb_database.py --rank-genes`) can.

### 5. Smaller differences from GTDB

- **Background genome.** It is random DNA at the species' GC content, ~0.15-0.4 Mb against GTDB's 3-5 Mb.
  Reads off the markers are fewer, and there are no paralogs or mobile elements to attract reads. A depth
  feature such as `sample_log_fragments` therefore sits on a different scale than at GTDB.
- **Gene order.** Without `--operons`, markers are shuffled per species and spaced evenly, so there are no
  neighbouring genes for pairs or long reads to span. The default `build_mini_db.sh` world and the tuning world
  have none.
- **Marker loss.** It is uniform across genes and clades. In GTDB some markers are missing in whole clades.

## Recommendations

1. **Fix the categoriser** (`simulate_gtdb_release.py`, `marker_category`): ribosomal for `u[SL]\d+`,
   `arch_S`, `L3_arch`, `rpl4p`; translation for `eif`, `aif`, `ef-1`, `erf1`, `pelota`, `rna_pol`,
   `kow_elon`; conserved for `sec[EY]` and `atp`. This changes only the `--gene_rates categories` worlds; the
   rates use their own generator, so the rest of a release stays the same.
2. **Better: take the rates from the real factors.** For example `--gene_rates r226`, reading each gene's
   within factor from a table shipped next to `markers_r226.tsv` (made from `gene_congeners.tsv`, as
   `rates_vs_r226.tsv` here). The mini worlds would then have GTDB's gene-to-gene pattern exactly, including
   archaea's, and no name heuristics.
3. **Real gene lengths.** Ship a mean length per marker (from the r226 `reference.map`) and use it instead of
   the HMM length.
4. **Treat archaeal conservation-feature results** from the 0.7.x benchmark world with care until 1 or 2 is
   done; bacterial results stand.
5. **Choose reduced gene sets only from real builds**, as `build_gtdb_database.py --rank-genes` already does.
   Do not tune the ranking on mini worlds.

Side observation: `local/v10/model_logs/relatives_by_gene_conservation.txt` traced 0 genes at r226 v10 (Spearman
"nan"): the trace took a read's genome from the mini worlds' contig names; fixed the same day
([2026-10-05-trace-relatives.md](../2026-10-05-trace-relatives.md)).
