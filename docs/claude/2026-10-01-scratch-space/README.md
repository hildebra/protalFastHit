# Disk space and file churn of the training data collection

2026-10-01, branch `audit-fixes` at `c575aa5` with the uncommitted `--scratch` and progress changes
to `scripts/build_gtdb_database.py` and `scripts/collect_training_data.py`. Question: a GTDB r226
build on the HPC sat for a long time with one `python3` (the collector) in state D at 0.7% CPU.
How much local scratch would the simulated samples need, if they were made on a node's own disk
(`--scratch`) instead of the network file system that holds OUTDIR?

## Measured

`measure.sh` (WSL, `simulate_metagenomes` built from the same tree, pbsim3 from
`envs/protal-db-build.yaml`, genomes of the synthetic 900-species world
`~/protal-perf/world900`, whose sequences compress like real ones; 50,000 pairs per setup,
pbsim3 at depth 20 of one 282 kb genome of 3 contigs):

| Reads | gz (kept) | plain (temporary) |
|---|---:|---:|
| 100 bp, HS20 | 194 bytes per pair | 471 bytes per pair |
| 150 bp, HSXt | 177 bytes per pair | 671 bytes per pair |
| 250 bp, MSv3 | 384 bytes per pair | 1071 bytes per pair |
| PacBio (`errhmm:ERRHMM-SEQUEL`) | 0.35 bytes per base (merged, gzip -1) | pbsim3: fq.gz 0.32, maf.gz 0.62 per base, `.ref` = the genome |
| Nanopore (`qshmm:QSHMM-ONT-HQ`) | 0.91 bytes per base | fq.gz 0.86, maf.gz 0.53 per base, `.ref` = the genome |

pbsim3 writes these three files for each sequence of the genome (`r_0001.*`, `r_0002.*`, ...).
SAMs: 93 MB of `.sam.zst` per million pairs of the gene-dense `w900` world
([load-and-output](../2026-09-30-load-and-output/README.md)); real genomes put a smaller share of
their reads on marker genes, so this is an upper bound.

## What the default design writes

From the code: the simulator writes a sample's reads plain (one FASTQ pair per genome from ART,
appended to the sample's plain R1/R2), with a decompressed copy of each of the sample's genomes,
then BGZF-compresses the sample and deletes the rest; the design points run in parallel
(`--threads` at a time, at most all 15 training and 18 test points), their samples one after the
other. For long reads the collector writes a plain copy of every genome of every sample of a
design point, pbsim3 adds its files, and the collector merges each sample's reads into one fq.gz
before it deletes the point's temporary folder; one point at a time.

Assumed: 4 MB per plain genome; ~150 genomes per training sample (20-200 species, 0.3 and 0.1
for a second and third strain) and ~260 per test sample (10-300, 0.5 and 0.2); about 100 contigs
per genome (drafts and MAGs), for the file counts only.

| | Kept until the end | Temporary, at its peak |
|---|---:|---:|
| training pe reads (22.5M pairs) | 5.7 GB | 12 GB (15 points at once: genome copies 9 GB, the 500k-pair samples' plain reads 2.2 GB) |
| test pe reads (15.2M pairs) | 3.8 GB | 24 GB (18 points: genome copies 19 GB, 1M-pair samples 4.4 GB) |
| training pb, ont reads (2.25 Gb each) | 0.8 + 2.1 GB | 18.5 GB (deepest ONT point: 1,800 genome copies twice 14.4 GB, pbsim3 and merged reads 4.1 GB) |
| test pb, ont reads (1.12 Gb each) | 0.4 + 1.0 GB | 10.6 GB |
| SAMs (pe, se, pb, ont; upper bound) | 2-7 GB | |
| profiles, dumps, tables | ~2 GB | |
| **sum** | **~18-23 GB** | |

The largest moment is the test set's paired-end simulation, with the training samples kept:
about 13 GB + 24 GB = **~40 GB**, up to ~55 GB when large communities and the deepest samples
coincide. With fewer threads than design points, the genome copies scale down. **60 GB of local
scratch is the least; 100 GB leaves room** for larger communities or a deeper design.

Temporary writes, deleted again: about 470 GB over the run (genome copies 410 GB, plain reads
55 GB), in roughly 140,000 files for the paired-end samples and, with ~100 contigs per genome,
several million for the long reads (each genome copy plus three pbsim3 files per contig). A
network file system is slow at creating and deleting that many files; a node's own disk is not.
That fits the collector seen in state D: in the long-read stage the Python process itself copies
the genomes, merges pbsim3's files and deletes them.

## Changed

`build_gtdb_database.py --scratch DIR` makes the training and test samples in `DIR/training` and
`DIR/test` (the simulators' temporary files are inside them) and copies the tables to
`OUTDIR/training` and `OUTDIR/test`; the trainer and the parity check use the samples there. The
run reports the space it takes on DIR every `--progress-every` seconds, after each collection, and
the most it took at the end (from `statvfs` every 5 s, less what DIR held at the start), so that
the next run can be sized from a real one.

Not done: the genome copies could be written once per run instead of once per sample (the
simulator decompresses every genome of every sample, the collector copies every genome of every
long-read sample), which would remove most of the 470 GB.

## Follow-up: the builds' `full_reference.fna`

Seen on the same HPC run: `protal_db/full_reference.fna` 86 GB next to a 21 GB
`database.protal`, and the training database has a near copy (without the species it leaves out).
The converter writes it (the marker genes of all genomes of `genomic_files_all`), `protal --build`
reads it for the uniqueness check and the gene conservation factors, and nothing reads it after
the build: a rebuild converts the release again, since the build packed `reference.fna`.

Changed: `build_gtdb_database.py` removes each database's `full_reference.fna` once its build is
done (and on a rerun that finds the build done), which frees ~170 GB of OUTDIR at r226. The
converter deletes each of its chunks as it joins them into the file, so the conversion no longer
holds the 86 GB twice.

Left in OUTDIR, not on scratch: the builds read the file front to back, which a network file
system does well, and making the conversion and both builds on scratch would need about 250 GB
more there (two ~86 GB files and two builds' outputs). Not measured: how much smaller a
zstd-compressed `full_reference.fna` would be (`--build` reads `.zst` inputs); a species' genomes'
copies of a gene are near-identical but not next to each other in the file.

## Follow-up: `full_reference.fna` compressed

`zstd_full_reference.sh` (WSL; zstd 1.5 CLI, 4 compression threads, 1 decompression thread) on the
GTDB-like tuning world's full reference (`~/tune/V3/protal_db/full_reference.fna`, 285 MB, 257,718
records, simulated strains):

| zstd | ratio | compress | decompress |
|---|---:|---:|---:|
| `-3` | 11.0 | 641 MB/s | 740 MB/s |
| `-3 --long=27` | 11.5 | 375 MB/s | 719 MB/s |
| `-6 --long=27` | 12.8 | 243 MB/s | 805 MB/s |
| `-9 --long=27` | 13.2 | 190 MB/s | 790 MB/s |
| `-12 --long=27` | 13.1 | 89 MB/s | 708 MB/s |
| `-19 --long=27` | 16.2 | 2 MB/s | 937 MB/s |

Changed: the converter writes `full_reference.fna.zst` with `zstd -6 --long=27 -T<threads, up to 8>`
(piped, so no raw copy is written first), and plain `full_reference.fna` only without the `zstd`
command; `--from_db` reads either and writes the copy the same way; `protal --build` reads the
`.zst` as it is (`zstd::InputFile`, also when given the raw name). At the world's ratio the 86 GB
of r226 would be ~7 GB, written in ~3 minutes on 8 threads, and each pass of a build over it
decompresses in ~2 minutes; GTDB's real ratio is not measured (the converter's log gives the
size).

## Follow-up: measured on the r226 build (2026-10-02)

SLURM job 23865669 (scripts `ff57266`, 16 threads, `--scratch` on the node's SSD with 157.4 GB
free; logs in `local/protal0.7_r226_v1/`, see
[2026-10-02-r226-build-evaluation](../2026-10-02-r226-build-evaluation/README.md)): the run took at
most **23.7 GB** on scratch (the estimate above: ~40 GB, up to ~55), and the samples left at the end
took 15.5 GB (estimate ~18-23 GB): 9.6 GB of training and 5.8 GB of test samples. GTDB r226's
`full_reference.fna.zst` was **5.32 GB for 79,520,648 sequences** (16x smaller than the 86 GB of
the plain file in the run before), the training database's copy 5.0 GB; both were removed after
their builds. The job read 2.4 TB and wrote 0.99 TB in all (sacct), most of it the simulator
reading every genome of the table at each design point (fixed in `c40cec3`).
