# protal 0.7 against the shipped 0.6.0a: speed, memory, detection, abundance, strains, database builds

- **Date**: 2026-09-30.
- **Code**: 0.7 = branch `audit-fixes` at `4b21427` (Version 0.7.0), built in WSL from a Linux-FS copy
  (`~/fix-build`, Release). 0.6.0a = the tag `0.6.0a` (`014f4a9`), exported with `git archive` and
  built the same way (`~/protal-0.6.0a`; it too says `v0.6.0`, as every build before the bump did).
- **Machine**: WSL2 Ubuntu 24.04 (kernel 6.18), Intel Core Ultra 7 258V, 6 threads given to WSL,
  23 GB, gcc 13.3. Quiet (load below 1) except for these runs, which ran one at a time with `-t 6`.
- **Data**:
  - The GTDB-like tuning world `~/tune/release_p` (765 species, 3 genomes each; the
    [tuning study](../2026-09-29-model-training-tuning/README.md)), converted once
    (`gtdb_to_protal_db.py`) and built by each version from the same files.
  - The independent test set of the [read-type models run](../2026-09-30-read-type-models.md)
    (`~/tune/V2/test`): 10–300 species per sample, abundances lognormal with sigma 2, second and third
    strains with probabilities 0.5 and 0.2, from the world's representative and other genomes. Its
    PacBio and Nanopore samples replay the paired-end communities (pbsim3).
  - One deep sample simulated here: 5M pairs, 2x150, 150 species, strains 0.3/0.1 (`deep.sh`).
  - The strain audit's world and run A (`~/audit5`, 8 species, 42 samples with one strain per species
    at 1–50x; [strain audit](../2026-09-29-strain-audit/README.md)), with its true genotypes.
- **Models**: 0.6.0a uses the model it ships (`random_forest.xml`, as `model.xml`). 0.7 uses the
  four models trained on this world by `build_gtdb_database.py` (`~/tune/V2/trained_model*.xml`,
  added with `--add_model`). `v07m06` is 0.7 with 0.6's model, to tell the binary's changes from the
  model's. The strain world's databases both use the shipped model.
- **Run**: `scripts/build.sh`, `profile.sh`, `diagnose.sh`, `heldout.sh`, `margin_sweep.sh`,
  `deep.sh`, `build_scaling.sh`, then `score.py` (detection, abundance, time, memory) and
  `eval_strains.sh` (the strain audit's `evaluate.py`, copied here with two changes: each version's
  qcmsa is re-run, and 0.6's 0-based partitions are read with an offset). Outputs in `~/bench07`;
  the tables are in `results/`.

## Summary

On the same data, 0.7 runs about twice as fast as 0.6.0a (5M pairs in 49 s instead of 100 s; 3.2×
faster on the strain run with qcmsa) in the same memory (3.4–3.6 GB, which the index sets), builds a
database 3× faster, and stores it in 107 MB instead of 3.4 GB. With models trained for the database it
detects species better, most at low depth (F1 0.74 against 0.14 at 500 read pairs), and its strain MSAs
recover true SNPs at low depth that 0.6 lost (SNP recall at covered true SNPs at 2x 0.65 → 1.00).

One regression: 0.7's abundances are further from the truth than 0.6's (Bray-Curtis 0.061 against
0.035 on the deep sample). The cause is `--depth_identity_margin 0.04`, which leaves out a species'
reads whose identity is more than 0.04 below that of its best reads (about 1.0): a 150 bp read with more
than 6 differences. Strains of this world differ from the reference by 0.4–4% at the markers, and single
reads scatter around that, so 0.04 drops 6–26% of the bases of strains 1–4.5% away (2x150; more for
100 and 250 bp reads). Species simulated from other strains are underestimated, and reference-genome
species overestimated in relative terms. With 0.08, 0.7's
abundances are the best of all, with every species in the database (0.016 on the deep sample) and
with whole clades and 20% of species missing from it; 0.08 still resists relatives' reads slightly
better than no margin. Recommendation: make 0.08 the default (checked on the round-3 world first).

## Speed and memory

Whole runs (index load, alignment, profiling, strain MSAs without qcmsa but for run A), mean of the
two samples of a point; wall time varies by up to 30% between identical runs at these sizes.

| Scenario | 0.6.0a wall | 0.7 wall | 0.6.0a CPU | 0.7 CPU | 0.6.0a peak RSS | 0.7 peak RSS |
|---|---|---|---|---|---|---|
| 5M pairs 2x150, 150 species | 100 s | 49.1 s | 310 s | 187 s | 3.40 GB | 3.56 GB |
| 500k pairs 2x250 | 22.9 s | 13.1 s | 66.0 s | 47.9 s | 3.39 GB | 3.56 GB |
| 500k pairs 2x150 | 10.2 s | 6.0 s | 30.4 s | 23.5 s | 3.39 GB | 3.58 GB |
| 500k pairs 2x100 | 10.0 s | 5.0 s | 23.7 s | 17.1 s | 3.39 GB | 3.57 GB |
| 10k pairs 2x150 | 2.4 s | 0.9 s | 3.0 s | 3.4 s | 3.30 GB | 3.46 GB |
| strain run A: 42 samples, 8 species, with qcmsa | 116 s | 36.6 s | 332 s | 144 s | 3.20 GB | 3.29 GB |

On the deep sample the index loads in 1.7 s instead of 9.4 s and the alignment takes 33.7 s instead of
74.3 s. 0.6 writes plain SAM; 0.7 writes `.sam.zst` in the alignment threads, within those times.

## Detection and abundance

F1 of the species reported (knob 0.5) and Bray-Curtis dissimilarity between the true and the reported
relative abundances (0 = identical; the simulator's `relative_abundance` summed per species; a
species missing on one side counts as 0). `m0.08` = 0.7 with `--depth_identity_margin 0.08` (the SAMs
of 0.7 profiled again; detection is unchanged by the margin, within 0.001).

| Scenario | 0.6.0a | 0.7 | 0.7 with 0.6's model | 0.7 m0.08 |
|---|---|---|---|---|
| 5M pairs 2x150 | 0.993 / 0.035 | 0.993 / 0.061 | – | 0.993 / 0.016 |
| 500k 2x250 | 0.979 / 0.034 | 0.996 / 0.070 | 0.987 / 0.073 | 0.995 / 0.015 |
| 500k 2x150 | 0.980 / 0.049 | 0.994 / 0.062 | 0.982 / 0.066 | 0.993 / 0.023 |
| 500k 2x100 | 0.984 / 0.059 | 0.993 / 0.093 | 0.977 / 0.121 | 0.993 / 0.029 |
| 10k 2x150 | 0.944 / 0.303 | 1.000 / 0.056 | 0.935 / 0.112 | – |
| 500 2x150 | 0.143 / 0.669 | 0.743 / 0.262 | 0.197 / 0.603 | – |

Species the database lacks (the database without the ~20% of species and whole clades that the
training database leaves out; scored on the species it has, their true abundances renormalised):

| Scenario | 0.6.0a | 0.7 | 0.7 m0.08 |
|---|---|---|---|
| 500k 2x250 | 0.972 / 0.044 | 0.985 / 0.067 | 0.985 / 0.020 |
| 500k 2x150 | 0.982 / 0.049 | 0.998 / 0.064 | 0.998 / 0.027 |
| 500k 2x100 | 0.983 / 0.070 | 0.998 / 0.099 | 0.997 / 0.037 |
| 10k 2x150 | 0.921 / 0.318 | 0.984 / 0.074 | 0.984 / 0.028 (margin 1) |

- The detection gain is the model's: with 0.6's model, 0.7 detects about as 0.6 does (0.977–0.987
  against 0.979–0.984 at 500k). 0.7's models were trained on samples of this world (another design and
  seed than the test set), 0.6's on older databases, so the gain is what a database trained for its
  reference gives, not a like-for-like comparison of models.
- `--whole_read_alignment` gives the same abundances as the default (Bray-Curtis equal to 4 digits):
  the bias is not the anchored alignment. `--depth_identity_margin 1` gives what 0.08 gives with all
  species in the database, and slightly worse with species missing (0.040 against 0.037 at 2x100).

The margin swept on the deep points (Bray-Curtis, mean of the two samples):

| margin | 0.04 | 0.06 | 0.08 | 0.10 | 0.15 | 0.25 |
|---|---|---|---|---|---|---|
| 2x100, all species | 0.093 | 0.041 | 0.029 | 0.029 | 0.029 | 0.029 |
| 2x100, species missing | 0.099 | 0.048 | 0.037 | 0.038 | 0.040 | 0.040 |
| 2x150, all species | 0.062 | 0.027 | 0.023 | 0.023 | 0.023 | 0.023 |
| 2x150, species missing | 0.064 | 0.028 | 0.027 | 0.028 | 0.028 | 0.028 |
| 2x250, all species | 0.070 | 0.024 | 0.015 | 0.015 | 0.015 | 0.015 |
| 2x250, species missing | 0.067 | 0.027 | 0.020 | 0.020 | 0.021 | 0.021 |

By kind of species (500k 2x250, sample 1): of the true abundance of species simulated from their
representative genome alone (0.180), 0.6 reports 0.183 and 0.7 0.222; of species with a representative
and another strain (0.329), 0.326 and 0.291. On 2x100 the species simulated from another strain alone
(0.369) get 0.366 and 0.316.

### What the margin keeps, by where the reads come from

The margin is in read identity as `AlignmentIdentity` (`Strain.h`) computes it: M / (M + X + I + D)
over the read's CIGAR, soft clips left out, so an indel's bases count as differences. The threshold is
the 98th percentile of the taxon's read identities, weighted by aligned bases (`TopIdentity`), minus the
margin; the 98th percentile is about 1.0 wherever some of a species' reads match the reference exactly.
0.04 therefore drops a 100 bp read with more than 4 differences, a 150 bp read with more than 6, and a
250 bp read with more than 10.

In `simulate_gtdb_release.py` every genome of a species, the representative included, is mutated from
the species' ancestral sequence with its own divergence (substitutions, a share of them codon indels),
drawn per genome from 0.002–0.02 for this world (`~/tune/world/simulation/divergence.tsv`). A strain
therefore differs from the representative, the database's reference, by the sum of two draws: 0.4–4%,
about 1 − ANI at the markers (in the simulation the markers diverge as fast as the rest of the genome).
A read's differences are roughly Poisson around its length times that divergence plus the sequencing
error rate: at 3.5% on 150 bp, 5.25 expected differences, and P(more than 6) = 0.28.

`scripts/identity_margin.py` applies protal's rule to 0.7's SAMs of the 500k-pair samples (sample 1,
primary alignments), and tells a taxon's reads apart by the genome in their simulated names. Share of
aligned bases kept, full database:

| Reads of | 2x100 at 0.04 / 0.08 | 2x150 at 0.04 / 0.08 | 2x250 at 0.04 / 0.08 |
|---|---|---|---|
| the reference genome itself | 0.995 / 1.000 | 1.000 / 1.000 | 0.959 / 0.999 |
| a strain 0–1% from the reference | 0.963 / 1.000 | 0.990 / 1.000 | 0.882 / 0.998 |
| a strain 1–2% | 0.879 / 0.998 | 0.935 / 1.000 | 0.785 / 0.998 |
| a strain 2–3% | 0.768 / 0.994 | 0.856 / 0.999 | 0.647 / 0.994 |
| a strain 3–4.5% | 0.636 / 0.981 | 0.744 / 0.995 | 0.598 / 0.992 |
| other species, on this taxon | 0.474 / 0.908 | 0.371 / 0.864 | 0.370 / 0.884 |

With species missing from the database, the other species' share of the taxa's aligned bases grows
from 1–4% to 5–20%, and 0.04 keeps 33–50% of them, 0.08 keeps 85–89%; the own-read shares are as above.
At 0.15 every read is kept. So the margin separates a species' strains from its relatives poorly at any
value in this world, whose congeneric species are 3–12% apart at the markers; 0.08 wins because a
species' own reads far outnumber its relatives'. 250 bp reads lose more at 0.04 than 150 bp reads,
even from the reference itself: ART's 250 bp profile (MSv3) has more errors towards the reads' ends.
Real strains of a species can differ from its representative by up to ~5% genome-wide (the 95% ANI
species boundary), so 0.04 would undercount them too.

### Follow-up: dynamic thresholds

A fixed margin below the 98th percentile cannot follow the strain a sample holds: its reads centre at
about 1 − (the strain's difference from the reference) − (sequencing errors). Alternative anchors, tried
with an experimental build (`depth_rule_patch.py`: the rule from `$PROTAL_DEPTH_RULE`; not in the
repository) on 0.7's SAMs of the deep points (`dynamic_margin.sh`). The gene median is the median over a
taxon's genes with 3 or more reads of each gene's median read identity: a strain covers all of a
species' genes, while relatives' reads pile on the conserved ones. `top:0.04` reproduces 0.7's numbers.

| Rule | 2x100 / 2x150 / 2x250 | same, species missing | 5M pairs |
|---|---|---|---|
| 98th percentile − 0.04 (0.7) | 0.093 / 0.062 / 0.070 | 0.099 / 0.064 / 0.067 | 0.061 |
| 98th percentile − 0.08 | 0.029 / 0.023 / 0.015 | 0.037 / 0.027 / 0.020 | 0.016 |
| gene median − 0.04 | 0.031 / 0.029 / 0.016 | 0.042 / 0.031 / 0.021 | 0.018 |
| gene median − 0.06 | 0.029 / 0.023 / 0.015 | 0.039 / 0.026 / 0.020 | 0.016 |
| gene median − 0.08 | 0.029 / 0.023 / 0.015 | 0.040 / 0.027 / 0.020 | 0.016 |
| gene median − 3 SD of a read's identity (binomial, at its divergence and length) | 0.036 / 0.040 / 0.018 | 0.046 / 0.040 / 0.023 | 0.023 |
| gene median − 5 SD | 0.028 / 0.027 / 0.014 | 0.039 / 0.029 / 0.020 | 0.019 |
| 98th percentile − (0.04 + 1 × (1 − gene median)) | 0.031 / 0.028 / 0.015 | 0.042 / 0.031 / 0.021 | 0.017 |
| 98th percentile − (0.04 + 2 × (1 − gene median)) | 0.028 / 0.024 / 0.015 | 0.039 / 0.027 / 0.020 | 0.016 |

(Bray-Curtis, mean of the two samples; detection unchanged.) By kind of species (`group_bias.py`), the
good rules remove the bias of 0.04 alike: on 2x100 sample 1, species from another strain alone (true
0.369) get 0.362–0.364, mixtures with the representative (0.437) 0.433–0.437, the representative alone
(0.110) 0.119–0.122; 0.04 gives 0.316, 0.458 and 0.152. A minor strain cut off by a gene median that
follows the dominant one did not show: this world's strains are at most 4% from the reference.

Anchoring on the gene median takes most of the error away even at 0.04, and from 0.06 on it equals a
fixed 0.08 below the top, as the rule that scales the margin by the gene median's divergence does. None
beats the fixed 0.08 here, because the world does not reach either end where they would differ: strains
near the species boundary (~5% from the representative, where 0.08 below the top starts to drop 100 bp
reads) and close strains beside relatives (where the gene median sits higher and cuts closer). That
0.08 and no margin at all score nearly alike with species missing (0.037 and 0.040 on 2x100) shows the
other reason: protal's depth is already a median over genes, which relatives' reads on a few conserved
genes hardly move.

## New read types (0.7 only)

0.6.0a profiles paired-end reads only. The same communities as single-end, PacBio and Nanopore reads:

| Scenario | F1 | Bray-Curtis | wall | peak RSS |
|---|---|---|---|---|
| single-end 500k 150 bp (first mates) | 0.994 | 0.063 | 4.3 s | 3.52 GB |
| single-end 10k 150 bp | 1.000 | 0.047 | 1.4 s | 3.46 GB |
| PacBio 90 Mb | 1.000 | 0.093 | 7.9 s | 3.50 GB |
| Nanopore 90 Mb | 1.000 | 0.090 | 10.3 s | 3.54 GB |
| Nanopore 3 Mb | 0.982 | 0.124 | 1.6 s | 3.49 GB |

## Strains

Run A, strain rows (a genome other than the reference), after qcmsa, by depth: the share of the
genomes' gene positions called, SNP recall at covered true SNPs, N at covered true SNPs, and false
alternative calls per million called invariant sites.

| depth | called 0.6 / 0.7 | SNP recall 0.6 / 0.7 | N at true SNPs 0.6 / 0.7 | false alt ppm 0.6 / 0.7 |
|---|---|---|---|---|
| 1x | 0.18 / 0.57 | 0.57 / 0.993 | 42.8% / 0.52% | 126 / 499 |
| 2x | 0.42 / 0.80 | 0.65 / 0.996 | 35.2% / 0.36% | 0 / 254 |
| 3x | 0.59 / 0.92 | 0.71 / 0.998 | 28.7% / 0.21% | 0 / 140 |
| 5x | 0.76 / 0.94 | 0.83 / 0.999 | 17.3% / 0.09% | 0 / 37 |
| 10x | 0.85 / 0.94 | 0.95 / 1.000 | 4.8% / 0.01% | 0 / 1 |
| 20x | 0.86 / 0.94 | 0.99 / 1.000 | 0.65% / 0% | 0 / 0 |
| 50x | 0.78 / 0.96 | 1.000 / 1.000 | 0.04% / 0% | 0 / 0 |

- 0.6 leaves many true SNPs as N at low depth (its strand filter and allele rules, fixed in rounds 1
  and 4), and its qcmsa removes more cells the deeper the sample (0.78 called at 50x, 0.92 before it). 0.7 writes positions with one read (step O), which costs
  false alternative calls at 1–3x.
- 0.7's called share at 10–20x is 0.94 over all strain rows but 0.99 for bacteria: its MSAs hold only
  the samples whose profile reports the species (`--msa_knob`), and the shipped model misses the
  archaeon in some samples (round 4). 0.6 writes the archaeon's rows regardless.
- 0.6.0a writes 0-based partition files (fixed in round 2), which puts every gene boundary one column
  off in IQ-TREE's partitions and in the evaluation unless corrected.

## Database builds

The tuning world (765 species, `reference.fna` 95 MB, `full_reference.fna` 285 MB), from the same
converted files; the conversion took 2.0 s at 6 threads, `--add_model` 0.7–1.5 s per model.

| Build | threads | wall | CPU | peak RSS | on disk |
|---|---|---|---|---|---|
| 0.6.0a, raw files | 6 | 40.3 s | 78.5 s | 3.47 GB | 3.82 GB |
| 0.6.0a, raw files | 1 | 31.0 s | 30.9 s | 3.47 GB | 3.82 GB |
| 0.7, raw files (`--no_compress`) | 6 | 13.1 s | 27.9 s | 3.65 GB | 3.80 GB |
| 0.7, raw files | 1 | 13.7 s | 13.7 s | 3.48 GB | 3.80 GB |
| 0.7, `database.protal` at zstd level 3 | 6 | 4.9 s | 21.1 s | 4.35 GB | 129 MB + full reference |
| 0.7, `database.protal` at level 19 (default) | 6 | 101.2 s | 180.5 s | 4.88 GB | 107 MB + full reference |

- 0.6.0a is slower on 6 threads than on one. Its uniqueness check was slower on 8 threads than on 4
  ([index build report](../2026-09-29-index-build-parallel.md)); 0.6 has no stage timers to show it here.
- At this size a raw build is bound by writing the 3.4 GB index; the single file at level 3 does not
  write it and takes 5 s. Level 19 spends 29 s compressing the index and 69 s packing the rest, for a
  file 17% smaller than level 3's.
- The strain world (8 species): 0.6.0a 6.6 s, 0.7 raw 2.3 s, 0.7 single file 8.0 s.
- At GTDB r226 the index build is estimated at 30–60 min on 16 cores
  ([index build gains](../2026-09-30-index-build-gains/README.md)); not measured here.

## Not covered

Real (non-simulated) samples, a GTDB-scale database, and the time of long-read alignment on reads
over 65 kb. The samples come from one synthetic world whose strains differ from the references by
0.4–4% at the markers; how the margin behaves on real strain divergences at the GTDB markers is the
thing to check before changing the default.
