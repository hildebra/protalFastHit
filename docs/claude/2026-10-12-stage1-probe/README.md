# What would a sparse first stage keep? The --stage1_probe option

- **Date**: 2026-10-12.
- **Question**: a two-stage protal would look the reads up in a sparse index first and align only the reads that
  stage keeps, against the taxa it found ([the proposal](../2026-10-12-profiler-fork-ideas.md), with its follow-up on
  keeping the candidate reads in a temporary file). How few reads can the first stage keep without losing reads
  that align, and which taxa lose reads? The same counts also answer the [seeding profile](../2026-10-12-seeding-profile/README.md)'s
  open question: can a mate whose first lookups match only by chance be given up?
- **Code**: the working tree on `dfbee7d` (with `e04f4a1`'s whole-block prefetch); not committed yet. The probe is a
  development option: a run with `--stage1_probe` writes the same SAM and profile as one without, and three more
  tables.
- **Status**: implemented and checked locally ([below](#local-check-the-archaea-world)); the r226 cluster run is
  pending (it needs the probe committed and pulled on the cluster).

## What a first-stage rule is

protal looks a read's k-mers up by their 15-base core. The 16 bases around it (the flex cell) are scored by the bases
they share with the read's, and a lookup's best score is that of its best cell (16: all of them). The first stage
would keep the cores whose hash falls in 1/c of its range (an index about 1/c the size, with every entry of those
cores), and keep a read if at least m of its lookups there score T or more. A pair counts both mates' lookups.

Because the sparse index has the same entries for its cores as the full one, the full index gives the first stage's
scores: the probe takes every lookup's best score as the seeding computes it (lookups the seeding skips once a read has
its seeds are scored apart) and its core's level, the trailing zero bits of a hash of the core (in the index of 1/c
for c = 2^level and every smaller c; so the sparser indexes are subsets of the denser). A key of one value has no flex
cells (a stage-1 index would have to store one): its score is the flank bases the read shares with the gene beside the
entry's core, in the orientation in which the read's core is the gene's.

At the read's end, after the output handler, each rule is tallied:

- c = 1, 2, 4, 8, 16, 32; T = 10 to 16; m = 1, 2, 3, 5;
- for all reads (pairs for paired-end reads) and for the reads protal writes records of: the share a rule keeps,
  and its **recall** (the share of the reads with records it keeps);
- per taxon of the best record (each gene's best for long reads): its reads with records and those each rule keeps
  (m = 1 and 2).

And per mate of a short read: the best score among its first k lookups in read order (k = 1, 2, 4, 8, all), by
whether the mate aligns, only its fragment does, or neither.

## The outputs

In the run's misc folder, per sample:

| file | rows |
|---|---|
| `<sample>_stage1_probe.tsv` | every rule: `c`, `min_score`, `min_hits`, the reads, those kept, `share_kept`, the reads with records, those kept, `recall` |
| `<sample>_stage1_probe_taxa.tsv` | per taxon of the best records: `taxid`, `name` (its lineage, as the profile writes it), the reads with records, then one column per rule (`c8_t14_m1`, ...) |
| `<sample>_first_lookups_probe.tsv` | short reads: per k and best score, the mates whose fragment does not align, whose other mate does, and that align |

The log has one line per sample with a few rules. [probe_summary.py](scripts/probe_summary.py) turns the tables and
the profile into tables of share kept and recall by c and T, the cheapest rule at a recall of 99%, 99.5% and 99.9%,
the called species that lose half or more of their reads under those rules, and the give-up shares.

## Running it at r226

[probe_cluster.sh](scripts/probe_cluster.sh) builds protal at the checkout's HEAD and runs the real paired-end sample
(B10, 50.6M pairs) and the HiFi barcode with `--stage1_probe` on 32 threads, then the summaries. It needs the probe
committed and pulled on the cluster:

```bash
cd ~/dev/protalFastHit && git pull
D=docs/claude/2026-10-12-stage1-probe/scripts
sbatch --exclusive -c 32 --mem 100G --wrap "bash $D/probe_cluster.sh $PWD \
    /hpc-home/hildebra/DB/protal/protal0.8.1_r226_v24/protal_db/database.protal \
    /hpc-home/hildebra/data/raw/Illumina/PEARL_AGE/00_fastq/B10_R1_001.fastq.gz \
    /hpc-home/hildebra/data/raw/Illumina/PEARL_AGE/00_fastq/B10_R2_001.fastq.gz \
    /hpc-home/hildebra/data/raw/PacBio/Data_Package_B20230626/Falk_Hildebrand_EI_FH_ENQ-5441_A_01_B20230626-Demultiplex_Barcodes/pb_demux_ccs_r64319e_20230622_123341-1_A01/outputs/m64319e_230622_130841.bc1002--bc1002.hifi_reads.fastq.gz \
    runs/stage1_probe"
```

## Local check: the archaea world

WSL, cores 0-3, the working tree built Release; the [seeding profile](../2026-10-12-seeding-profile/README.md)'s local
world (every r226 archaeal species' marker genes, `~/seed/db`) and its first 100,000 simulated pairs of 150 bp (reads
of whole genomes of 120 species with strains and congeners, so most come from outside the marker genes); 4 threads,
`--no_strains`. Tables in [results/local_archaea](results/local_archaea/).

- **Nothing changes**: with and without `--stage1_probe`, the SAM records (as a sorted set) and the profiles are the
  same, for the pairs, for mate 1 as single-end reads, and for mate 1 as PacBio reads (`--no_profile`). ctest
  525 (3 new: `Stage1Probe.*`), all pass.
- **A first version counted keys of one value as score 16**: in this world many cores have one value, so 99.98% of
  the pairs had a "perfect" lookup at c = 1. Such lookups are now scored against the gene (above).
- **The numbers** (not r226's: this world has 43% of the cores in its index and 4 cells a block, r226 99% and 70):
  2.89% of the pairs have records. T = 14, m = 1: c = 1 keeps 3.02% of the pairs at a recall of 100%, c = 2 2.96%
  at 99.20%, c = 4 2.85% at 97.09%, c = 8 2.68% at 92.03%, c = 16 2.31% at 79.74%
  ([pe_summary.md](results/local_archaea/pe_summary.md)). Single-end: c = 4 2.40% at 96.53%. So here the full
  index's near-exact lookups separate the reads well (3% kept, nothing lost), and the sparse indexes lose reads
  faster than the estimate of 44/c own-taxon hits per pair predicted (c = 8: 8% lost, not 0.4%).
- **Giving up mates** (first 4 lookups at most 11): 97.5% of the mates would be given up, and 6.0% of the mates
  that align and 19.4% of the mates of aligned fragments would be lost; over all of a mate's lookups, 0.7% of the
  aligned mates have none above 11 (they align through their other mate).

**How to read it**: the cheapest rule at a recall of 99.5% or more gives F, the share of pairs a first stage would
write to its temporary file, and c the size of its index (about 26.6 / c GB at r226). The called species that lose half
their reads or more are what a two-stage protal would risk; low-depth species are expected among them. For giving up
mates: a threshold that loses no aligned mate (or only mates whose fragment aligns by the other) changes no call.
