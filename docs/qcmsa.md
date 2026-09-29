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

It writes (the prefix defaults to `<species>`, the input path without `.raw.msa.fna`):

| Output | What it is |
|---|---|
| `<prefix>.msa.fna` | the filtered MSA, the one to use downstream |
| `<prefix>.partition.txt` | the partition with recomputed coordinates |
| `<prefix>.qcmsa_summary.tsv` | machine-readable decision log: what was removed and why (`--no-summary` skips it) |
| `<prefix>.qc.png` | optional MRate2 heatmap (`--plot`, needs matplotlib) |

## What it filters, in order

1. Coverage gate: drops genes and gap-fills (sample, gene) cells that are too sparsely covered,
   from the meta columns `hcov` and `mean_vcov_nonzero`. protal itself writes every observed gene
   and sample, so this is the only place these thresholds exist.
2. Multi-allelicity (MRate2) filter: removes genes and samples that are multi-allelicity outliers
   by an iterative Tukey-IQR rule, and masks single outlier cells.
3. Site cleanup: drops uninformative variable sites (and, if asked, constant sites).
4. Sequence floor (`--reapply-hcov`): drops whole sequences with too few valid bases.

## Parameters

### Multi-allelicity (contamination, mixed strains)

| Flag | Default | Effect |
|---|---|---|
| `--preset strict\|default\|sensitive` | default | sets the two below: strict = 1.0 / 1, default = 1.5 / 2, sensitive = 2.0 / 3 |
| `--iqr-mult FLOAT` | 1.5 | Tukey fence multiplier; lower removes more |
| `--min-bad INT` | 2 | bad peers a gene or sample needs before removal; higher removes less |
| `--sample-abs-min-bad INT` | 0 (off) | remove a sample that is multi-allelic in at least this many genes, whatever the fence says; 1-2 catches mixed or conspecific strains |
| `--gene-abs-min-bad INT` | 0 (off) | remove a gene that is multi-allelic in at least this many samples |
| `--mrate2-include-zeros` | off | compute the fence over all items, zeros included, so that on a clean baseline any multi-allelic item is flagged |
| `--max-mrate2 FLOAT` | from the data | hard per-cell MRate2 cap for masking outlier cells |
| `--no-mask-cell-outliers` | masking on | do not mask single outlier cells |

Explicit `--iqr-mult` or `--min-bad` override the preset.

### Coverage

| Flag | Default | Effect |
|---|---|---|
| `--gene-min-hcov FLOAT` | 0.3 | minimum fraction of a gene covered for a cell to pass; 0 disables |
| `--gene-min-mean-depth FLOAT` | 1.0 | minimum mean depth over covered positions; 0 disables |
| `--gene-min-samples INT` | 3 | drop a gene unless more than this many samples pass coverage; 0 disables. Only samples in the MSA count |

### Sites and sequences

| Flag | Default | Effect |
|---|---|---|
| `--min-parsimony-samples INT` | 2 | drop variable sites where fewer than N samples differ from the majority; 0 keeps all variable sites |
| `--discard-constant` | off | drop constant sites (kept by default: they inform branch lengths) |
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
      --preset sensitive --gene-min-hcov 0.3 --gene-min-mean-depth 1 --gene-min-samples 3
```

The `.qcmsa_summary.tsv` says what was removed (`genes_filtered_coverage`,
`genes_filtered_mrate2`, `samples_filtered`, `coverage_gap_filled_cells`, `sites_in` to
`sites_kept`, and the reasons per gene and sample).

For all species of a run at once, `just strain-refilter` (see the [justfile](../justfile) and
[testing.md](testing.md#strain-test-harness)) loops over a run's raw MSAs:

```bash
just strain_variant=test2 strain-refilter refilter_hcov=0.3 refilter_depth=1 refilter_min_samples=3 preset=sensitive
```

## Picking values

- Too much removed, strains collapse: loosen with `--preset sensitive` (higher `--iqr-mult`, higher
  `--min-bad`), lower or disable the coverage gate, and set `--min-parsimony-samples 0`.
- Contamination or mixed strains slip through: tighten with `--preset strict`, or add
  `--sample-abs-min-bad 1`.
- Your own filtering downstream: run protal with `--no_qcmsa` and filter the `.raw.msa.fna` with
  the `.meta.tsv`, which carry everything qcmsa uses.

The SNP filters (`--snp_*`) cannot move to qcmsa: they need per-read base qualities and strands,
which the meta table does not have, so they stay in protal.
