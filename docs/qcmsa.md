# qcmsa: filtering the strain MSAs

`qcmsa` (`scripts/qcmsa.py`, installed as `qcmsa`) post-filters protal's per-species strain MSAs,
using the per-(sample, gene) metrics protal writes to `.meta.tsv`. protal runs it by default after
writing each raw MSA, but it is also a standalone tool, so you can re-filter with other parameters
without re-running protal. The [website](https://protal.earlham.ac.uk/main.php?site=documentation#strain-resolution)
introduces the strain outputs; this page is the operator's guide.

## Inputs and outputs

For one species, qcmsa reads three files from the strain output directory:

| Input | What it is |
|---|---|
| `<species>.raw.msa.fna` | protal's own MSA, the input to filter |
| `<species>.raw.partition.txt` | gene to column ranges (RAxML style, 1-based) |
| `<species>.meta.tsv` | per-(sample, gene) coverage and multi-allelicity metrics, the filter's evidence |

The first row of protal's MSA, `<species>_reference`, is the database's genome of the species.
qcmsa keeps it and does not treat it as a sample: it is not in `.meta.tsv`, the sequence floor
does not apply to it, and it does not count towards `--min-parsimony-samples`. IQ-TREE sees it as
one more leaf; in the 2026-09-29 strain audit it neither helped nor hurt the trees.

It writes (the prefix defaults to `<species>`, the input path without `.raw.msa.fna`):

| Output | What it is |
|---|---|
| `<prefix>.msa.fna` | the filtered MSA, the one to use downstream |
| `<prefix>.partition.txt` | the partition with recomputed coordinates |
| `<prefix>.qcmsa_summary.tsv` | machine-readable decision log: what was removed and why (`--no-summary` skips it); when no MSA is left, its `status` line says why |
| `<prefix>.qc.png` | optional MRate2 heatmap (`--plot`, needs matplotlib) |

qcmsa first removes these outputs of an earlier run, so none of them outlives a run that writes it
no more. It stops with an error if the MSA names a sequence twice.

## What it filters, in order

1. Coverage gate: drops genes and gap-fills (sample, gene) cells that are too sparsely covered:
   the share of the gene that the sample's row writes (from the MSA itself), and optionally the
   gene's mean depth (from `.meta.tsv`). protal itself writes every observed gene and sample, so
   this is the only place these thresholds exist.
2. Multi-allelicity (MRate2) filter: removes samples, then genes, whose multi-allelic rate is an
   outlier, and masks single outlier cells. See below.
3. Sequence floor (`--reapply-hcov`): drops whole sequences with too few valid bases in the genes
   left.
4. Site cleanup: drops columns without any A/C/G/T and, if asked, constant sites and
   low-parsimony sites. It judges A/C/G/T only: an IUPAC code is an ambiguity, as IQ-TREE reads it.

## Parameters

### Multi-allelicity (contamination, mixed strains)

| Flag | Default | Effect |
|---|---|---|
| `--preset strict\|default\|sensitive` | default | sets the two below: strict = 1.0 / 1, default = 1.5 / 2, sensitive = 2.0 / 3 |
| `--iqr-mult FLOAT` | 1.5 | Tukey fence multiplier; lower removes more |
| `--min-bad INT` | 2 | multi-allelic genes a sample needs (samples a gene needs, positions a cell needs) before removal; higher removes less |
| `--mrate2-min-rate FLOAT` | 0.004 | floor under every fence: nothing with a rate at or below it is removed or masked |
| `--sample-abs-min-bad INT` | 0 (off) | remove a sample that is multi-allelic in at least this many genes, whatever the fence says |
| `--gene-abs-min-bad INT` | 0 (off) | remove a gene that is multi-allelic in at least this many samples |
| `--max-mrate2 FLOAT` | from the data | hard per-cell MRate2 cap for masking outlier cells |
| `--no-mask-cell-outliers` | masking on | do not mask single outlier cells |
| `--mrate2-include-zeros` | | no effect, kept for old command lines |

Explicit `--iqr-mult` or `--min-bad` override the preset.

A sample's rate is its multi-allelic positions (the IUPAC codes of the MSA) over its positions
with at least 2 reads, pooled over its genes (`multi_allelic` and `counts_vcov2` in `.meta.tsv`).
A sample is removed when its rate is above the Tukey upper fence of all samples' rates
(Q3 + `--iqr-mult` × IQR, needing 4 samples) and above `--mrate2-min-rate`, and it is
multi-allelic in `--min-bad` genes. Genes are then judged the same way on the samples kept, and
cells (one sample, one gene, `multi_rate_vcov2`) on the samples and genes kept. There is one pass.
A pooled rate does not grow with depth, as a count of multi-allelic genes does; in the 2026-09-29
strain audit, two-strain mixtures had rates of 0.45–1.9% and single strains at most 0.34%, before
`--snp_min_af` and the read identity margin lowered the latter. Without the floor, a fence over
mostly clean samples is 0 and removes any sample with two IUPAC codes.

### Coverage

| Flag | Default | Effect |
|---|---|---|
| `--gene-min-hcov FLOAT` | 0.3 | minimum share of a gene's positions (its columns where the reference row has a base) that the sample's row writes as a base or IUPAC code, not `-` or `N`, for the cell to pass; 0 disables |
| `--gene-min-mean-depth FLOAT` | 0 (off) | minimum mean depth over the whole gene, from the reads the MSA takes (`hcov` × `mean_vcov_nonzero` in `.meta.tsv`) |
| `--gene-min-samples INT` | 1 | drop a gene unless more than this many samples pass coverage, so by default a gene needs 2 (as protal needs 2 samples for an MSA); 0 disables. Only samples in the MSA count |

### Sites and sequences

| Flag | Default | Effect |
|---|---|---|
| `--min-parsimony-samples INT` | 0 (off) | drop variable sites where fewer than N samples (not counting the reference row) differ from the majority. A site where one sample differs is that strain's own mutation: with 2, terminal branches shrink to 0–9% of their length (2026-09-29 strain audit), so use it for topology only |
| `--discard-constant` | off | drop constant sites (at most one of A/C/G/T). A tree then needs an ascertainment correction (IQ-TREE `+ASC`, e.g. `just strain-trees strain_iqtree_model=GTR+G+ASC`); without it, branch lengths come out tens of times too long |
| `--reapply-hcov INT` | 0 (off) | drop sequences with fewer than N valid (non-gap, non-N) bases; protal passes its `--msa_min_hcov` (default 1000) |
| `--partition-base auto\|0\|1` | auto | coordinate base of the input partition; auto detects it (protal before this version wrote 0-based files) |

## Setting qcmsa options from a protal run

protal passes only `--prefix` and `--reapply-hcov` (from `--msa_min_hcov`) itself; everything else
keeps qcmsa's defaults. Every qcmsa flag is reachable through `--qcmsa_args`, forwarded verbatim:

```bash
protal --map map.tsv --db DB -t 16 --qcmsa_args "--preset strict --gene-min-hcov 0.5"
protal --map map.tsv --db DB -t 16 --no_qcmsa          # skip qcmsa: only the raw MSAs
```

The arguments are appended last, so they also override what protal passes, e.g.
`--qcmsa_args "--reapply-hcov 0"` beats the value from `--msa_min_hcov`. protal does not check
them. If qcmsa fails for a species (a bad flag, say), or cannot be found at all, protal reports
it, keeps the raw MSA and exits with status 1 at the end of the run. If qcmsa keeps no gene or
sample, protal warns that there is no filtered MSA.

## Re-filtering an existing run

Run qcmsa on the files protal already wrote. Nothing is re-aligned; this takes seconds.

```bash
qcmsa out/strains/s__Bacteroides_ovatus.raw.msa.fna \
      out/strains/s__Bacteroides_ovatus.raw.partition.txt \
      out/strains/s__Bacteroides_ovatus.meta.tsv \
      --prefix out/refiltered/s__Bacteroides_ovatus \
      --preset sensitive --gene-min-hcov 0.3 --gene-min-mean-depth 1 --gene-min-samples 2
```

The `.qcmsa_summary.tsv` says what was removed (`genes_filtered_coverage`,
`genes_filtered_mrate2`, `samples_filtered`, `coverage_gap_filled_cells`, `sites_in` to
`sites_kept`, and the reasons per gene and sample).

For all species of a run at once, `just strain-refilter` (see the [justfile](../justfile) and
[testing.md](testing.md#strain-test-harness)) loops over a run's raw MSAs:

```bash
just strain_variant=test2 strain-refilter refilter_hcov=0.3 refilter_depth=1 refilter_min_samples=2 preset=sensitive
```

## Building trees

protal lists the species of each run in `<strain output dir>/species.tsv`: species, taxid, the
samples admitted, and the file names of the raw and filtered MSAs (`-` for none). Outputs of other
species in the directory are an earlier run's.
`just strain-trees` builds an IQ-TREE tree for each species in that list (`strain_tree_input=raw`
uses the raw MSAs), with `GTR+G` (`strain_iqtree_model`), 1000 ultrafast bootstraps and a fixed
seed (`strain_iqtree_seed`). It skips MSAs with fewer than 4 sequences, the reference row included.
Partitioned trees (`-p <prefix>.partition.txt`) came out the same as unpartitioned ones in the
2026-09-29 strain audit, at several times the run time.

## Picking values

- Too much removed, strains collapse: loosen with `--preset sensitive` (higher `--iqr-mult`, higher
  `--min-bad`), or lower or disable the coverage gate.
- Contamination or mixed strains slip through: tighten with `--preset strict`, or lower
  `--mrate2-min-rate`.
- Your own filtering downstream: run protal with `--no_qcmsa` and filter the `.raw.msa.fna` with
  the `.meta.tsv`, which carry everything qcmsa uses.

The SNP filters (`--snp_*`) cannot move to qcmsa: they need per-read base qualities and strands,
which the meta table does not have, so they stay in protal.
