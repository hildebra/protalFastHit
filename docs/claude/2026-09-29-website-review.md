# Website documentation review

2026-09-29 · Claude Code, for Falk Hildebrand

protal.earlham.ac.uk (Home, Documentation, Downloads, FAQ, as served on 2026-09-29) compared with
branch `audit-fixes` at `e917af3` (protal 0.6.0). Most of the website still holds. The items below
are wrong or missing for the current code; most changed on `audit-fixes`, so they matter once that
branch is released. Method: the pages were downloaded and read in full; claims were checked
against the source, the `--help`, `--full_help` and `--map_help` of a fresh WSL build of `e917af3`,
and a run of `examples/mini_db/run.sh` with that build (all checks passed).

"Changed in" names the commit or audit step (see [the code audit](2026-09-26-code-audit.md)) where
known.

## Documentation page

| Section | The website says | Now | Changed in |
|---|---|---|---|
| Useful options; Strain resolution | `--strain_preset strict\|default\|sensitive` | the option does not exist; use `--qcmsa_args "--preset strict"` (or `sensitive`) | `alpha` `014f4a9` |
| Useful options | `--knob`: no range given | must be 0 to 1, checked at start; choose it on data like yours (the old help's claim that 0.4-0.6 hardly changes F1 did not hold) | step H `d35f158` |
| Useful options | (missing) | `--depth_identity_margin` (default 0.04): only reads within this identity of a species' best reads count towards its abundance; reads of relatives still count for detection | step I `f6b7ba6` |
| Useful options | `--profile_only`: profile existing SAMs | outputs go to `-o` if given, else next to each SAM | step K `0d04d2c` |
| Useful options | `--model`: alternative PMML file | also takes the name of a model in the database folder; protal checks the model before aligning (exit 2 if unusable) | step H |
| Mapping file | `SAMPLEID`: "protal does not actually read this column"; names come from `PREFIX` | `#SAMPLEID` names the sample in MSA rows, `meta.tsv`, statistics and logs | step K |
| Mapping file | give each sample its own SAM and PROFILE, "otherwise the samples overwrite each other's output" | protal stops before starting if two samples share a file. Every row needs a value in every column of the header; an empty cell is an error | steps B `969b241`, K |
| Run a single sample | (not said) | with `-1/-2 -o DIR`, the SAM and profile files go directly into `DIR`; `strains/` and `misc/` are subfolders. A single sample gets no strain MSAs (see below) | |
| Strain resolution | MSAs are written for each species | only for species that pass the model in at least 2 samples (or those named by `--msa_species`), with rows only for samples where the species passes | step I |
| Output: alignments | `x.sam.gz`: "writes an uncompressed sam and then compresses it with pigz" | only when the SAM name ends in `.gz`; the example maps name `.sam` files and get plain SAMs. Either change the maps to `.sam.gz` or say so. SAMs are written as `.partial` and renamed when complete, so a crashed run is not reused | step 3 `ced6a2c` |
| Output: alignments | `x.sam.gz.err`: "diagnostic output" | `<SAM name>.err`: the reads whose alignment does not fit the database (gene missing, past a gene's end, bases differ); protal warns with their number | step K |
| Output: profiles | `x.profile.log`: no header, columns listed in prose | a header row: `Predicted Probability RepGenome Lineage Abundance VCovStdDev GeneVariance GeneVariance5 Name TaxID Summary VCov LowIdentityShare MeanGeneCov MeanGeneCovRatio GeneCov0 ...`; abundance is 0 for rejected taxa | steps I, K |
| Output: profiles | `x.profile.gene.log`: per-gene detail | header: `Sample TaxID Lineage Name GeneID Gene Reads CoverageSum GeneLength ANISum MAPQSum MeanANI MeanMAPQ CoveredBases CoveredFraction Mono Bi Tri Tetra FilteredNoAllele FilteredMono FilteredBi FilteredTri FilteredTetra` | step K |
| Output: profiles | `x.profile.genes.log` header `Truth, Predicted, ..., TaxaxAbundance, ..., UniqueTwoMersReads, UniqueTwoMerReads, ...` | `Predicted Probability TaxID Lineage TaxVCOV TaxAbundance GeneID GeneRefLength TotalReads TotalMappedLength MAPQ UniqueMers UniqueTwoMers UniqueMerReads UniqueTwoMerReads ANI VCov VCovExp HCovExp HCovObs HCovObsRel Consistency`; MAPQ and ANI are means, not sums | step K |
| Output: profiles | `x.profile.truth_annotated`: TP/FP annotation of the profile | the model's training dump: one row per species with reads, with the truth, the call, the probability and every feature. Also written from a map's `PROFILE_TRUTH` column | steps D `050e6cd`, H |
| Output: profiles | example `head` shows `d__Bacteria\|p__...` and `RS_GCF_...` | the current code joins the lineage with `;`; re-run the example and paste fresh lines | |
| Output: strains | `y.partition.txt`, `y.raw.partition.txt` | now 1-based, so `iqtree2 -p` accepts them (0-based files failed with "Negative site ID") | step C `ca2df11` |
| Output: strains | `y.snp_stats.tsv` columns | also a `refs_retained` column | step 1 `2da2748` |
| Output: strains | (missing) | `y.qc.png`, qcmsa's MRate2 heatmap, with `--qcmsa_args "--plot"` | |
| Output: misc | `hcov.tsv`, `snps_*.tsv`: "the remaining 120 columns" per gene | one column per marker gene id of the database: 168 for GTDB r226 (bac120 and ar53 together) | |
| Output: misc | `y.statistics.tsv` | written for single-sample species too; header `#SAMPLEID Accepted VerticalCoverage TotalReads TotalLength MeanAni MeanMAPQ` | step K |
| Output: misc | (missing) | `<prefix>_seedsizes_histogram.tsv`, `<prefix>_anchorsizes_histogram.tsv`, `<prefix>_runtime.tsv`: seeding and alignment diagnostics | |
| Output: misc, SNP filters | quality thresholds `--snp_min_phred_sum 90`, `--snp_min_mean_qual 15` | unchanged defaults, but base qualities are now read as Phred+33 (they were 3 too low), so slightly more variants pass. Deletions are called; N is never an allele | steps 1, C |
| Usage | (missing) | exit status: 0 = all done; 1 = finished, but a sample or output failed (listed at the end); other = stopped before aligning (bad options, inputs, database or model). Useful for Snakemake and Nextflow users | steps 3, A `5c93647` |
| Usage | "we will use data from .... (context)" | placeholder; name the study the four SRR samples come from | |
| Usage | "The mapping for for all these samples", "acccession" | typos | |
| Installation: compile from source | "Instructions are provided in the github Readme" | the README is now short; link `docs/installation.md` on GitHub instead | this change |
| Species profiles | `protal_profile_utils` "is installed alongside protal" | `conda-recipe/build.sh` and `just install` in this repository install `protal_map_utils` and `qcmsa` but not `protal_profile_utils`. Check the bioconda recipe; install it there, or tell users to take it from `scripts/` | |
| Download the database | the database extracts to a folder of files | databases built with this version are one file, `database.protal`; `--db` takes it or its folder. The r226 0.5.1a folder still loads. When a new database is published, update the text and the "compatible with" notes on the downloads page | `66016e9` |

## FAQ

| Question | The website says | Suggested change |
|---|---|---|
| Can I use a custom database? | no; GTDB r214 | yes: a database can be built from any GTDB release, or part of one, with `scripts/build_gtdb_database.py` (`docs/building-a-database.md`). The shipped database is r226; update the release and its species count |
| Is there a windows version? | no | still no native version, but it builds and runs under WSL2 (checked on Ubuntu 24.04 for this review) |
| Can I use protal on long reads? | short reads only | also: paired-end only; single-end reads are not supported yet |
| Can I track strains across samples? | "strain-trees and pairwise sample distances for each species" | protal writes no pairwise distance matrix at present (the similarity matrix is switched off in `StrainWrapper2`, `src/RunProtal.h`); distances come from trees built on the MSAs |
| What if my sample contains conspecific strains? | protal detects and discards them; "minimum and maximum distance to other samples" | the removal is qcmsa's multi-allelicity filter (`--qcmsa_args "--preset strict"` or `"--sample-abs-min-bad 1"`); the distances are not a protal output (see above) |
| Is protal using alignments or k-mers? | "WFA2 library (include link)"; "o excel" | link https://github.com/smarco/WFA2-lib; typo |
| Which aligner; alignment approach | "a proprietary aligner" (twice) | protal is open source (GPL-2.0); "proprietary" reads as closed. Say "its own aligner" or "a purpose-built aligner" |
| Why not bwa-mem or bowtie? | "(link to preprint)" | link https://doi.org/10.64898/2026.08.03.742433 |

## Not on the website, now in `docs/`

Rather than adding these to the website, they are documented in the repository; the website could
link to them.

- Building from source, static builds, the AVX2 launcher, WSL2: `docs/installation.md`
- Output locations per mode, reruns, exit status, the options `--help` shows but the website does
  not, alignment and developer options, environment variables: `docs/running.md`
- The database files, `database.protal`, compression and conversion: `docs/database-files.md`
- Building a database from GTDB, reduced marker sets, build and train in one command:
  `docs/building-a-database.md`
- Training the presence model: `docs/model-training.md`
- qcmsa's parameters and re-filtering: `docs/qcmsa.md`
- `simulate_metagenomes`: `docs/simulation.md`

## Found in the repository while checking

Not website issues, but found on the way; not changed here.

- `protal --help` for `-o` says SAMs and profiles go to the subfolders `alignments/` and
  `profiles/`; in `-1/-2` mode they go directly into `-o` (`src/Options.h`, the prefix and SAM
  defaults after `OptionsFromArguments`).
- `--no_qcmsa`'s help says protal "fails gracefully with a warning" if qcmsa cannot be found; it
  reports an error and exits 1 (`RunQCMSA`, `src/RunProtal.h`).
- The removed `scripts/qcmsa.md` documented `--remove-constant` / `--keep-constant`; qcmsa's flag
  is `--discard-constant`. `docs/qcmsa.md` has the right flags.
- The CI workflow no longer builds the mini database or runs the end-to-end tests (removed in
  `e917af3`); the old README still said it did.
- `scripts/build_gtdb_database.py` trains on simulations of species that are all in the database,
  so the model sees none of the novel-species negatives that `docs/model-training.md` recommends.
