# Strain MSAs and trees

For each species found in two or more samples, protal writes a multiple sequence alignment (MSA) of
its marker genes with a row per sample, from which a strain phylogeny can be built. The
[website](https://protal.earlham.ac.uk/main.php?site=documentation#strain-resolution) introduces
the outputs. This page covers how they are made, how to filter them, and how to build trees:

1. protal calls each sample's variants per species and writes the raw MSA
   (`<species>.raw.msa.fna`, its partition and `.meta.tsv`);
2. `qcmsa` filters it (`<species>.msa.fna`), by default as part of the run;
3. you build trees from the filtered MSAs.

## Which samples and genes enter an MSA

- **Samples.** A species gets an MSA if it passes the model in at least two samples, with one row
  per such sample, so a run of one sample writes no MSAs. A sample passes with a probability of at
  least its knob (the knob of the profile: `--knob`, the model's knob for the sample's depth, or its
  calibrated calls with `--fdr`; never a species the singleton rule vetoes) and of at least
  `--msa_knob` (by default the same knob). Lower `--msa_knob` to add samples the profile leaves out.
  `--msa_species s__Genus_species,...` writes MSAs for the named species only.
- **Genes.** An MSA takes the genes with reads in its samples. For a species with relatives in the
  database, genes without any long unique k-mer are left out where the species has such k-mers on
  90% of its genes or more: a relative may share such a gene unchanged, and its reads would show as a
  second strain. Foreign genes (genes whose reads come from another genome,
  [running.md](running.md#alignment-options)) are left out unless `--keep_foreign_genes`.
- **Reads.** A read enters the rows only within `--msa_identity_margin` (0.04) of the species' best
  reads, stricter than for the abundance (0.08): a relative's reads within 0.08 put 10 times more
  false calls into the raw MSAs of mixed strains.
- **Bases.** A position is written with `--msa_min_depth` (1) reads showing one allele; a second
  allele (an IUPAC code) needs `--snp_min_cov` reads and the other SNP filters, and disagreeing reads
  otherwise give `N`. A sample's row stays only with `--msa_min_hcov` (1000) bases that are not `N`
  or gaps.

The first row, `<species>_reference`, is the database's genome of the species. `species.tsv` in the
strain output folder lists the species of the run: taxid, samples admitted, and the raw and filtered
MSA file names (`-` for none). Outputs of other species in the folder are an earlier run's.

## Strain MSAs over several runs

An MSA holds the samples of one protal run. To build MSAs over the samples of several runs (studies
aligned apart, samples added later), profile their SAM files together: no read is aligned again and
no index is loaded, and the MSAs are those of one run over all the reads (byte-identical in the
[2026-10-07 check](claude/2026-10-07-sam-combine/README.md)). The SAMs must come from the same
database; protal checks their `@SQ` genes against it.

```bash
protal --db DB -t 16 -o combined --profile_only 'studies/*/alignments/*.sam.zst'
protal --db DB -t 16 -o combined --profile_only study1/alignments/a.sam.zst,study2/alignments/b.sam.zst
```

- **Patterns.** An item with `*`, `?` or `[...]` stands for the `.sam`, `.sam.gz` and `.sam.zst`
  files it matches, sorted; a run's `.err`, profile and `.partial` files are passed over. Quoted,
  protal expands it and says how many SAMs it matched; unquoted, the shell does, which works as
  well. A pattern that matches no SAM, or a folder, stops protal before it starts.
- **Sample names.** A sample is named after its SAM file (`a.sam.zst` is `a`). Where two files share a
  name (a folder per sample, each with `aln.sam.zst`, or `sa.sam.zst` in two studies), every sample is
  named by the folders in which the paths differ: `studies/S1/aln.sam.zst` and
  `studies/S2/aln.sam.zst` are `S1` and `S2`, `p/study1/alignments/sa.sam.zst` and
  `p/study2/alignments/sa.sam.zst` are `study1_sa` and `study2_sa`; protal notes it. `--prefix`
  (one per SAM, in the order of the expanded list) names them otherwise.
- **Outputs.** Each sample's profile is written again, into `-o`, with the MSAs in `-o/strains`.
  Without `-o` the profiles go next to the SAMs and `strains/` and `misc/` into the current folder.
  protal stops before it starts where that would overwrite an earlier run's profiles or strain list
  (`species.tsv`); `--force` writes them again. The `<sam>.err` file is always written next to its SAM.
- **Maps.** `protal_map_utils merge --map run1.map run2.map --out combined > all.map` writes one
  map of the runs' samples, in which every SAM a run wrote is named by its absolute path: `protal --map
  all.map` profiles them where they are and aligns only the samples without one (`--new-sams` gives
  every sample a new SAM, to align all of them again). A map written by hand works the same way: an
  absolute path in its `SAM` column may point into any run's folder, and the read files of a sample
  whose SAM exists need not exist (protal warns).

Re-profiling costs a fraction of a run: at GTDB r226, profiling took 3.7 s of a 42 s paired-end run
([report](claude/2026-10-06-performance-profiling/README.md)). Until the MSAs are written, every sample
keeps what its MSA rows need of the species that enter MSAs, packed: per gene, its alleles at variant
positions (32 bytes each) and its coverage (a byte per base below 256x); the reads' records and the
other taxa are freed once its profile is written (`--taxon_statistics` keeps the taxa). A dense sample
of 1M read pairs from 60 species at about 25x keeps 37 MB (159 MB up to 0.7.8), so memory grows by that
per sample: 16 such samples on 4 threads peaked at 1.4 GB (2.8 GB up to 0.7.8). The MSAs are the same.

## Strains in long-read samples

A PacBio or Nanopore sample whose reads show two or more strains of a species gets a row per strain,
`<sample>_hap1`, `<sample>_hap2`, ... (the most abundant first), instead of one row with IUPAC codes
wherever both strains' alleles pass. A long read covers several marker genes and shows, at each
multi-allelic site, the allele of the one strain it comes from:

1. Sites that reads link, directly or through other sites, form blocks.
2. In each block, haplotypes grow from the reads: one starts from the read with the most sites and
   takes the reads that agree with it; a read that fits none starts the next.
3. Blocks are joined by their haplotypes' shares of the reads, which are the strains' shares of the
   sample in every block. A join needs to be at least 20 times likelier than any other.
4. Each row is called from its strain's reads with the same SNP filters, so it carries its strain's
   alleles also where the sample's filters did not see them, and `-` where its reads do not reach.

Phasing needs phased blocks on at least 3 genes: a second allele on one or two genes is another
genome's copy of a gene, not a strain. Reads are cut between two genes that are unlikely neighbours
(a chimera), and strains of about equal abundance, or too few reads, leave a block unphased.
`<species>.haplotypes.tsv` records each block (genes, sites, reads, the haplotype of each row and the
log odds of its join). `--no_phasing` writes one row per sample; paired-end and single-end samples
always get one ([report](claude/2026-10-02-phasing-and-foreign-genes/README.md)).

## Options

| Option | Default | |
|---|---|---|
| `--no_strains` | off | no MSAs or SNP tables. Variants are still called, since the model uses them, so profiles are the same |
| `--msa_species` | all | species to write MSAs for |
| `--msa_knob` | `--knob` | the probability a sample's species needs to enter the MSA |
| `--msa_min_hcov` | 1000 | bases (not `N`, not gaps) a sample's row needs to stay; passed to qcmsa as `--reapply-hcov` |
| `--msa_min_depth` | 1 | reads a position needs to be written, else `-` (2 before the 2026-09-29 strain audit) |
| `--msa_identity_margin` | 0.04 | how far below the species' best reads (98th percentile of identity) a read may be to enter the rows; 1 lets every read in |
| `--snp_max_alleles` | 3 | alleles an IUPAC code may encode: 1 only the top allele, 2 two-allele mixtures (R, Y, ...), 3 also three-allele ones (B, H, ...) |
| `--no_phasing` | off | one row per long-read sample |
| `--keep_foreign_genes` | off | keep foreign genes in depth and MSAs |
| `--no_qcmsa`, `--qcmsa_args`, `--qcmsa_script` | | skip qcmsa; extra qcmsa arguments ([below](#running-qcmsa-from-protal)); the qcmsa executable ([installation.md](installation.md#installing-a-source-build)) |

The SNP filters (`--snp_min_cov`, `--snp_min_phred_sum`, `--snp_min_mean_qual`, `--snp_min_af`,
`--snp_no_strand`) are described on the website. Base qualities are Phred+33. `--snp_min_af` is
0.15 (0.2 for Nanopore reads): an allele needs that share of a position's reads to be called or to
enter an IUPAC code. The profiles do not depend on these filters. They need per-read qualities and
strands, which `.meta.tsv` does not carry, so they cannot move to qcmsa.

## Filtering with qcmsa

`qcmsa` (`scripts/qcmsa.py`, installed as `qcmsa`) filters a species' raw MSA by the per-(sample,
gene) coverage and multi-allelicity that protal writes to `.meta.tsv`. protal runs it after each raw
MSA; it is also a standalone tool, to filter again with other values without running protal.

| Input (from the strain folder) | |
|---|---|
| `<species>.raw.msa.fna` | protal's MSA |
| `<species>.raw.partition.txt` | gene to column ranges (RAxML style, 1-based) |
| `<species>.meta.tsv` | per-(sample, gene) coverage and multi-allelicity |

| Output (`<prefix>`: the input without `.raw.msa.fna`, or `--prefix`) | |
|---|---|
| `<prefix>.msa.fna` | the filtered MSA, the one to use downstream |
| `<prefix>.partition.txt` | the partition with recomputed coordinates |
| `<prefix>.qcmsa_summary.tsv` | what was removed and why; when nothing is left, its `status` line says why (`--no-summary` skips it) |
| `<prefix>.qc.png` | with `--plot`: a heatmap of the multi-allelic rates (needs matplotlib) |

qcmsa removes an earlier run's outputs first, and stops if the MSA names a sequence twice. The
reference row is kept and is not a sample: it is not in `.meta.tsv`, the sequence floor does not
apply to it, and it does not count for `--min-parsimony-samples`.

### What it filters, in order

1. **Coverage.** Genes and (sample, gene) cells that are too sparsely covered: the share of the gene
   a sample's row writes, and optionally the gene's mean depth. protal writes every observed gene and
   sample, so these thresholds exist only here.
2. **Multi-allelicity** (MRate2): samples, then genes, whose rate of multi-allelic positions is an
   outlier, then single outlier cells are masked.
3. **Sequence floor** (`--reapply-hcov`): sequences with too few valid bases in the genes left.
4. **Sites.** Columns without A/C/G/T, and optionally constant and low-parsimony sites. IUPAC codes
   count as ambiguities, as IQ-TREE reads them.

A sample's rate is its multi-allelic positions over its positions with 2 or more reads, pooled over
its genes. It is removed when the rate is above the Tukey fence of all samples' rates (Q3 +
`--iqr-mult` × IQR, from 4 samples on), above `--mrate2-min-rate`, and it is multi-allelic in
`--min-bad` genes; genes are then judged the same way on the samples kept, and cells on both. A
pooled rate does not grow with depth. In the 2026-09-29 strain audit's simulations single-strain rows
stayed at 0.09% or less, while mixtures with a minor strain of 15% or more reached 0.34% or more. The
fence finds outliers: where most samples of a species are mixtures it rises above them, and only
`--sample-abs-min-bad` removes them.

### Parameters

| Flag | Default | Effect |
|---|---|---|
| `--preset strict\|default\|sensitive` | default | `--iqr-mult` / `--min-bad`: 1.0 / 1, 1.5 / 2, 2.0 / 3. Explicit values override it |
| `--iqr-mult` | 1.5 | Tukey fence multiplier; lower removes more |
| `--min-bad` | 2 | multi-allelic genes a sample needs (samples a gene needs, positions a cell needs) before removal |
| `--mrate2-min-rate` | 0.002 | nothing at or below this rate is removed or masked |
| `--sample-abs-min-bad`, `--gene-abs-min-bad` | 0 (off) | remove a sample (gene) multi-allelic in at least this many genes (samples), whatever the fence |
| `--max-mrate2` | from the data | per-cell cap for masking |
| `--no-mask-cell-outliers` | | do not mask single cells |
| `--gene-min-hcov` | 0.3 | share of a gene's positions a sample's row must write (a base or IUPAC code) for the cell to pass; 0 disables |
| `--gene-min-mean-depth` | 0 (off) | minimum mean depth over the gene (`hcov` × `mean_vcov_nonzero`) |
| `--gene-min-samples` | 1 | a gene needs more than this many samples passing coverage (so 2 by default) |
| `--reapply-hcov` | 0 (off) | drop sequences with fewer valid bases; protal passes `--msa_min_hcov` |
| `--min-parsimony-samples` | 0 (off) | drop sites where fewer samples differ from the majority. With 2, terminal branches shrink to 0-9% of their length, so use it for topology only |
| `--discard-constant` | off | drop constant sites. Trees then need an ascertainment correction (IQ-TREE `+ASC`), else branch lengths come out tens of times too long |
| `--partition-base auto\|0\|1` | auto | coordinate base of the input partition (protal before 0.7 wrote 0-based ones) |

What to change:
- **Too much removed, strains collapse**: `--preset sensitive`, or loosen the coverage gate.
- **Contamination or mixed strains slip through**: `--preset strict`, or a lower `--mrate2-min-rate`.
- **Your own filtering**: run protal with `--no_qcmsa` and use the `.raw.msa.fna` and `.meta.tsv`,
  which carry everything qcmsa uses.

### Running qcmsa from protal

protal passes `--prefix` and `--reapply-hcov` (from `--msa_min_hcov`); `--qcmsa_args` appends any
other flags, last, so they also override those two:

```bash
protal --map map.tsv --db DB -t 16 --qcmsa_args "--preset strict --gene-min-hcov 0.5"
protal --map map.tsv --db DB -t 16 --no_qcmsa          # only the raw MSAs
```

protal does not check the arguments. If qcmsa fails for a species or cannot be found, protal reports
it, keeps the raw MSA and exits with status 1 at the end. If qcmsa keeps nothing, protal warns.

### Filtering an existing run again

Nothing is aligned again; it takes seconds:

```bash
qcmsa out/strains/s__Bacteroides_ovatus.raw.msa.fna \
      out/strains/s__Bacteroides_ovatus.raw.partition.txt \
      out/strains/s__Bacteroides_ovatus.meta.tsv \
      --prefix out/refiltered/s__Bacteroides_ovatus \
      --preset sensitive --gene-min-hcov 0.3 --gene-min-mean-depth 1 --gene-min-samples 2
```

`just strain-refilter` does it for every species of a run
([development.md](development.md#strain-test-harness)).

## Building trees

`just strain-trees` builds an IQ-TREE tree per species of `species.tsv` (`strain_tree_input=raw` for
the raw MSAs): `GTR+G` (`strain_iqtree_model`), 1000 ultrafast bootstraps, a fixed seed
(`strain_iqtree_seed`). MSAs with fewer than 4 sequences, the reference row included, are skipped.
In the 2026-09-29 strain audit the reference row neither helped nor hurt the trees, and partitioned
trees (`-p <prefix>.partition.txt`) came out as unpartitioned ones at several times the run time.
