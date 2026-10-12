# The false positives and false negatives of r226 v24b

**Data.**
- The r226 v24b build ([2026-10-12-r226-v24b](../2026-10-12-r226-v24b/README.md); protal 0.8.2 `4a42f3b`, SLURM
  24171430, share archive unpacked to `local/v24b`): its models' calls (`model_logs/trained_model*.calls.tsv.gz`), the
  training tables (`work/{training,test}/training_data*.tsv`), and the reads behind every error
  (`model_logs/error_reads/<read type>/taxa.tsv.gz`, written by `scripts/error_reads.py` during the build). For
  every false positive and false negative this table lists where the taxon's reads came from
  (`fragments_on_taxon_from`). For every false negative it also lists where the taxon's own reads went:
  `own_best_on_taxon`, `own_best_elsewhere(_on)`, `own_unaligned`.
- The errors are at the knob protal calls at (the calls table's `call`: pe 0.72, se 0.7, pb and ont 0.5). Training
  samples are scored with species held out, and the test set by the final model. The v24b report's tables use knob
  0.5, so its counts differ a little.
- **Scripts** (WSL, cores 0-3, `~/micromamba/envs/protal-db-build`; outputs beside them):
  [error_anatomy.py](error_anatomy.py) → `error_anatomy_<type>.txt` (where the errors are, the false positives'
  sources, the false negatives by genome and fragments); [deep_dive.py](deep_dive.py) → `deep_dive_<type>.txt`
  (where missed taxa's reads went, the MAPQ filter, recurring pairs, per sample); [causes.py](causes.py) →
  [causes.txt](causes.txt) (one cause per error and its F1 cost); [twin_separability.py](twin_separability.py) →
  `twin_separability_{pe,pb}.txt` (taxa beside a held-out congener: features and AUC); [own_reads.py](own_reads.py) →
  [own_reads.txt](own_reads.txt) (the missed taxa's own reads).

**Question (the user's, 2026-10-12).** Analyse the false positives and false negatives of this run.

## Summary

- **The errors are in the scenario samples.** The 30 moderate, soil and soil_shallow samples (a thousand or more
  species each) hold 84% of pe's training false positives and 85% of its false negatives, and about 90% for long
  reads. The 309 design samples alone reach pe F1 0.9708 (all samples: 0.9560), PacBio 0.9836 (0.9714).
- **Most false negatives cannot be found from their reads.** 60% of pe's and 76% of the long reads' missed taxa have
  at most 2 reads of their own whose best alignment is on them; 13-34% have none. Such a species is present by the
  simulation's truth, but its marker genes received too few of its reads. A taxon has a row only because a
  congener's reads put one there. Counting as present only the taxa with 3 or more own reads would raise F1 by up to
  +0.012 to +0.017 (an upper bound: called taxa with fewer own reads would leave too).
- **60-80% of the false positives are a held-out congener's reads.** These are clouds of 3 or more fragments (26-40%),
  1-2 stray reads (20-40%), or twins within 0.02 on the marker genes (12-17%). Another genus's reads cause 4-12%,
  and a present congener's 3-8%. 31% of pe false positives come from 270 recurring taxa. The top pairs are GTDB twins
  0.011-0.03 apart, such as Priestia aryabhattai beside Priestia zanthoxyli, called in 12 of its 24 absent rows.
- **A present species beside a held-out congener is the false positives' mirror.** About 23% of the false negatives
  are "swamped": most reads on the taxon are its held-out congener's. Their own reads are a median 3 (pe) and 0
  (PacBio). They look *more* like a congener's cloud than the false positives do: identity AUC 0.17 (pe) and 0.19
  (pb), missed versus falsely called. No feature in the tables separates them, because the reads do not.
- **The errors are confident.** 44% of pe false positives score 0.9 or more, and 52% of the false negatives below 0.3.
  No knob trades one for the other.
- **What is left to win from the reads is small.** pe false negatives with 3 or more own reads are 40% of them. They
  are mostly the MAPQ filter's losses (16% of pe false negatives, at most +0.004 F1) and divergent strains below 0.97
  identity (9% pe, 18% ONT).

## 1. Where the errors are

Training samples (species held out), at protal's knob ([error_anatomy_*.txt](error_anatomy_pe.txt), section 1):

| | design samples | FP / FN in design samples | FP / FN in the scenario samples (pe, pb, ont 50; se 30) | F1 design | F1 all |
|---|---|---|---|---|---|
| pe | 309 | 406 / 502 | 2,431 / 3,342 | 0.9708 | 0.9560 |
| se | 309 | 418 / 556 | 1,161 / 1,975 | 0.9676 | 0.9610 |
| pb | 185 | 90 / 104 | 1,116 / 1,604 | 0.9836 | 0.9714 |
| ont | 185 | 137 / 130 | 1,453 / 1,928 | 0.9780 | 0.9649 |

Per sample (pe, medians): design samples have 41 present species, 2.0 false positives per 100 present and a 2.7%
miss rate; soil 2,019 present, 4.2 and 6.5%; soil_shallow 1,715 present, 6.3 and 6.9%; moderate 1,384 present, 2.5 and
3.4%. Among design samples the rates hardly change between 31 and over 100 present species. The scenarios' extra
errors come from what such communities hold: many congeners per genus, and in soil 60% of the species
missing from the training database.

## 2. False positives

One cause per false positive, from the main source of the reads on the taxon (at least half of them;
[causes.txt](causes.txt)), training samples:

| cause | pe | se | pb | ont | F1 cost, pe |
|---|---|---|---|---|---|
| cloud of a held-out congener (same genus, 3+ fragments, farther than 0.02) | 36% | 35% | 26% | 26% | 0.0064 |
| a held-out congener's stray reads (1-2 fragments) | 20% | 20% | 39% | 40% | 0.0035 |
| twin of a held-out species (within 0.02) | 13% | 12% | 17% | 15% | 0.0023 |
| a held-out species of another genus | 16% | 14% | 10% | 10% | 0.0028 |
| another genus's reads (a database species) | 10% | 12% | 4% | 5% | 0.0019 |
| a present congener's reads | 6% | 7% | 3% | 4% | 0.0011 |
| false positives | 2,837 | 1,579 | 1,206 | 1,590 | |

- **Fragments on the taxon** (pe): 1 for 35%, 2 for 13%, 3-10 for 27%, 11-100 for 20%, more for 5%. The clouds
  carry the reads (median identity 0.979 against 0.993 for the other sources): many reads at 2-3% divergence, the
  species boundary as the marker genes see it.
- **Other genera.** pe has 294 false positives, three quarters of them a single read, at identity 0.993
  from a database species of another genus. These are
  conserved or transferred genes, and the depth feature leaves them to the deep samples. Another 442 come from a
  held-out species of another genus. GTDB's suffixed genera (Bacillus_AD beside Priestia) are often close relatives.
- **Recurring taxa.** 270 pe taxa with 3 or more false positives hold 1,132 of 3,696 (training and test). 150
  (taxon, source) pairs with 3 or more hold 579 of the 2,837 training false positives. The top pairs are twins of
  GTDB r226 at marker distance 0.011-0.03: Pyrobaculum oguniense beside P. arsenaticum (9), Jeotgalicoccus
  nanhaiensis beside J. halotolerans (8), Cuniculiplasma sp021816385 beside sp027316895 (8), Priestia aryabhattai
  beside P. zanthoxyli and P. megaterium (12 false positives in 24 absent rows). In this build the source species is held out, so
  every sample that holds it puts its reads on the twin. In real use such a false positive appears only when the
  true species is missing from GTDB.
- **Scores.** 77% of pe false positives score 0.8 or more, 44% 0.9 or more, 25% 0.95 or more. For long reads it is half
  and a third. They are not near the knob.

## 3. False negatives

### What the reads hold

The missed taxa's own reads whose best alignment is on the taxon, any MAPQ ([own_reads.txt](own_reads.txt)):

| | missed | no own read | 1-2 own reads | 3 or more | F1 counting only 3+ as present (upper bound) |
|---|---|---|---|---|---|
| pe training | 3,844 | 18% | 42% | 40% | 0.9560 → 0.9707 |
| se training | 2,531 | 13% | 44% | 44% | 0.9610 → 0.9742 |
| pb training | 1,708 | 34% | 41% | 25% | 0.9714 → 0.9838 |
| ont training | 2,058 | 28% | 49% | 24% | 0.9649 → 0.9797 |
| pe test | 1,606 | 16% | 47% | 38% | 0.9557 → 0.9732 |

A taxon without own reads still has a row because its congener's reads (or others') landed on it. Its features are
those of an absent congener-holder, and the model rejects it as it should. Most of their reads were never on a
marker gene. Of the missed taxa's 2.3 M fragments in the SAMs, 94-96% did not align: reads of the rest of the genome.
The upper bound would be lower in a real recount, because called taxa with fewer than 3 own reads would leave the
truth too.

### Causes

One cause per false negative, the first that applies ([causes.txt](causes.txt)), training samples:

| cause | pe | se | pb | ont | F1 cost, pe |
|---|---|---|---|---|---|
| 1-2 fragments before and after the MAPQ filter | 36% | 38% | 42% | 45% | 0.0094 |
| swamped by a held-out congener (most reads on the taxon are its) | 24% | 22% | 22% | 23% | 0.0062 |
| the MAPQ filter (≤ 2 fragments after it, 3 or more before) | 16% | 14% | 8% | 8% | 0.0043 |
| divergent (identity below 0.97) | 9% | 11% | 9% | 18% | 0.0024 |
| own reads best on another species (half or more) | 4% | 4% | 3% | 3% | 0.0010 |
| swamped by another database species | 4% | 3% | 2% | 2% | 0.0009 |
| other | 8% | 7% | 14% | 1% | 0.0021 |
| false negatives | 3,844 | 2,531 | 1,708 | 2,058 | |

- **By genome** (pe): real strains 51% of the false negatives (1,975), in-silico strains 31%, representatives 17%. A
  present taxon with 1 fragment is missed 12% (representative) to 24% (real strain with alleles) of the time, with 3-10
  fragments 3-8%, with more than 100 fragments 0.6-2.4% ([error_anatomy_pe.txt](error_anatomy_pe.txt), section 3).
- **Representatives missed with 3 or more fragments** (281 in pe) are nearly all swamped. A median 12.5% of the reads
  on them are their own; the rest are a held-out twin's. Nanobdella sp028275775 had 250 own reads under 1,921 of
  Nanobdella stetteri_A's ([deep_dive_pe.txt](deep_dive_pe.txt), section 3).
- **Real strains lose reads to congeners.** Of a missed real strain's own aligned reads, a median 25% aligned best to
  another species (in-silico strains 5.5%, representatives 0%). Those reads went almost all to a congener of the same
  genus that was absent and not called (18,745 fragments in pe for strains with alleles). 1,268 went to an absent
  congener that was called: one strain's miss became another species' false positive.
- **The MAPQ filter**: 629 pe false negatives had 3 or more fragments before it and at most 2 after, 69 of them more
  than 10 before (their EM fragments a median 10.7). `em_fragments` is a feature already, and it gained nothing at
  r226 before (features.md).
- **Scores.** 21% of pe false negatives score below 0.05, 52% below 0.3, 71% below 0.5.

### Beside a held-out congener: the two errors are one population

Taxa with 3 or more fragments and a held-out congener within 0.15 in the sample, pe training
([twin_separability_pe.txt](twin_separability_pe.txt)): 6,416 called present taxa, 1,033 missed (800 swamped),
1,098 falsely called, 47,872 rightly rejected.

| feature (pe) | AUC present / absent | AUC missed / falsely called | median: swamped FN | FP | TP |
|---|---|---|---|---|---|
| identity | 0.947 | 0.174 | 0.957 | 0.975 | 0.992 |
| lu_per_kb | 0.906 | 0.261 | 24.7 | 45.4 | 63.5 |
| top_identity | 0.855 | 0.413 | 1.000 | 1.000 | 1.000 |
| excess_scaled_median | 0.063 | 0.791 | 0.034 | 0.014 | 0.000 |
| excess_high_share | 0.054 | 0.796 | 0.66 | 0.39 | 0.09 |
| low_identity_share | 0.446 | 0.670 | 0.072 | 0 | 0 |
| the model's p | 0.980 | 0.000 | 0.18 | 0.88 | 1.00 |

Every feature that tells present from absent points the wrong way between the two errors. The missed taxa look
more absent than the falsely called ones, and PacBio is the same (identity 0.188, excess_scaled_median 0.837). Their
own reads explain why: a swamped false negative has a median 3 own reads in pe (163 of 800 have 10 or more) and 0 in
PacBio (7 of 359). Only `low_identity_share`, 0.07 against 0 for the false positives, notices a second population of
reads. On 13 fragments that is one or two reads. The earlier error-read study ([r226 v15](../2026-10-07-error-read-signatures/README.md))
found no read signal at the species boundary; this repeats it for the case where both species are in the sample.

## 4. What this means

1. **The headline F1 is set by the scenario samples and by species the marker genes cannot see.** The design
   samples alone are at 0.971 (pe) to 0.984 (pb). Of the remaining false negatives, 60-76% are taxa without 3 reads of
   their own. Reporting F1 also with a truth that requires a few own reads (`error_reads.py` already counts them)
   would separate the method's errors from the data's limit. The build's console could print both.
2. **The false positives beside held-out congeners and the swamped false negatives are one population.** A feature
   or knob that lowers one raises the other: the v24b report's refits showed the training mix moving errors between
   the two at the same F1. What would help is not deciding the species: such a taxon could be reported as an unknown
   species of its genus. The r226 v15 study estimated +0.002-0.003 F1 for that, mostly in soil, with 82-92%
   of the added genus calls right in simulation (64% in gut;
   [2026-10-07-error-read-signatures](../2026-10-07-error-read-signatures/README.md)).
3. **The MAPQ filter's 16% (pe)** is the largest part that the reads do show. It is a ceiling of +0.004 F1, and earlier
   attempts (`em_fragments`) gained nothing. A soft count (reads below MAPQ 4 weighted by the EM's share) in
   `fragments` itself, instead of beside it, would be the next try. It is worth it only if the swamped taxa do not
   take the gain back as false positives.
4. **Divergent strains are the ONT-specific part** (18% of ONT's false negatives below identity 0.97, against 9% for
   PacBio). This and the borderline ONT rate loss in the v24b report point the same way: ONT alignment (Singletrack's
   ties, the indel bound) on divergent strains is worth a look.
5. **Recurring twins** (Priestia, Pyrobaculum, Jeotgalicoccus, Cuniculiplasma) are GTDB species pairs 1-3% apart on
   the markers. A list of such pairs (`species_clouds.tsv` already has them) could let protal report "X or its twin
   Y" when the evidence is the shared cloud. That trades a false positive for a less specific true call, as point 2
   does at genus level.
