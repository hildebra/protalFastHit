# The r226 v19 build: a rerun of v18, and how accurate the unknown share is

**Data.** The r226 v19 build (SLURM 24099010, q512n23, 2026-10-08 21:19 to 2026-10-09 01:09, 84 threads). Its protal,
simulator and scripts were at `a81cee3` (`build_metadata.tsv`), seed 1. Its share archive
`local/protal0.7.9_r226_v19_share.tar.gz` (516 MB) is unpacked to `local/v19`. It is compared with v18
([2026-10-08-r226-v18](../2026-10-08-r226-v18/README.md), `local/v18`). The composition numbers come from the build's
`model_logs/composition_accuracy*.tsv` (`scripts/composition_accuracy.py`, one line per sample), grouped with awk:

    cmp local/v18/model_logs/trained_model.varimp.tsv local/v19/protal0.7.9_r226_v19/model_logs/trained_model.varimp.tsv
    zcat .../trained_model.test_predictions.tsv.gz | md5sum      # v18 and v19
    # per design point, per depth bin and per share of the species the training database lacks, samples without a host:
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $c["host_share"]<0.01 {...}' model_logs/composition_accuracy_pb.tsv

**Questions (the user's, 2026-10-09).** Is there a noticeable degradation in performance, and where does it come from
(in theory for now)? How does protal predict the unknown fraction, and what influences it?

## Summary

- **v19 did not degrade: its models are v18's.** Every F1, FP and FN row of the summary equals v18's to the last digit.
  The variable importances are byte-identical, and so are the pe test predictions (same md5). v19 ran at `a81cee3`,
  which is before the strain alleles (`504de80`). Between v18 and v19 only these changed: the foreign scan became
  opt-in (it was in no model), the `--outdir` layout, the composition step, and the ancestry report's leave-one-genome-out.
  None of them reaches the profiles or the training tables. The build is also reproducible across nodes at 84 threads.
- **The step down is from v17, and it is the foreign leak** (v18 report). Against v15, the last honest build before
  it:
  - pe test F1 is level (0.9645 → 0.9622);
  - species held out is lower (0.9519 → 0.9462), and so are the soils (0.9655 → 0.9541; shallow soil 0.9396 → 0.9269)
    and gut (0.9962 → 0.9795);
  - the long reads gained (pb test 0.9785 → 0.9816, ont 0.9723 → 0.9794).
- **What explains the drop against v15, in theory:**
  - **Harder evaluation since v17:**
    - scenarios varied in richness and depth;
    - a moderate scenario;
    - strains 0.5,0.2, so more real strains, which are 67% of the FN;
    - in v18, a hold-out steered by species clouds, which is a new draw of held-out species.
  - **A lower knob:** 0.74 against 0.8 adds about 40 FP to the pe test set at about the same FN.
  - **Gut, unexplained:** it is the one drop not accounted for by a harder design (50 FN in 4 hold-out samples; v17's
    honest refit had 0.9871). It needs the error reads.
- **The ancestry report is honest now.** With a read's own genome left out of its species' alleles (`1c184fb`), the
  fixed sites' AUC falls from 0.81-0.86 to 0.58-0.68 (table in section 2). They still beat the congener sites by
  0.10-0.18. The polymorphic-site signal is real but a third to a half of what v18 reported.
- **The unknown share `?` is accurate for short reads from about 50,000 fragments up** (within ±0.02 of its truth).
  Below that it is too low:
  - −0.04 at 5,000-50,000 fragments;
  - −0.11 at 1,500-5,000;
  - −0.27 below 1,500.

  The called species of a shallow sample are the ones whose reads happened to land on their marker genes, so their
  depth is overestimated (+50% to +130%).
- **Long reads have a depth bias of their own, so their `?` is off at any depth:**
  - **PacBio:** depth +5% to +35%, `?` about −0.13 even in deep samples;
  - **Nanopore:** depth −12% to −19% in deep samples, `?` +0.04 to +0.10.

  Neither depends on the share of the sample's species that the database lacks. The PacBio cause is open.
- **A host's reads count as unknown genomes** of the called species' average size, which puts `?` at 0.85-0.92 in host
  samples. That is by definition, not an error of the estimate, but such a `?` is not a share of microbial cells.

## 1. v19 is v18 again

| | v18 | v19 |
|---|---|---|
| commit | `4b13640` → `e24dcff` (moved during the run) | `a81cee3` |
| feature set | `...+ancestry+gaps+untried` | the same |
| pe F1: species held out / test / soil / shallow soil / gut | 0.9462 / 0.9622 / 0.9541 / 0.9269 / 0.9795 | the same |
| pe test FP/FN | 106/243 | 106/243 |
| varimp, test predictions | | byte-identical |
| ancestry report: fixed-site AUC pe/se/pb/ont | 0.809 / 0.839 / 0.860 / 0.814 (circular) | 0.664 / 0.679 / 0.583 / 0.599 |
| composition check | not run | first run |

`504de80` (strain alleles) and `156ccda` were not on the HPC checkout. The strain alleles are still untested at r226.

## 2. The ancestry report without the circularity

The ancestry report's AUC, per taxon with 10 or more sites, FN taxa against FP taxa:

| | pe | se | pb | ont |
|---|---|---|---|---|
| identity | 0.513 | 0.547 | 0.387 | 0.471 |
| species' base at the congener sites | 0.562 | 0.571 | 0.407 | 0.461 |
| species' base at the fixed sites, v19 (v18) | 0.664 (0.809) | 0.679 (0.839) | 0.583 (0.860) | 0.599 (0.814) |
| fixed-site identity, v19 (v18) | 0.620 (0.779) | 0.652 (0.821) | 0.568 (0.831) | 0.556 (0.685) |

The own genome was left out in 19,064 FN records (pe), none of the FP's. The polymorphic-site mask separates missed
strains from novel congeners somewhat better than the congener sites alone (+0.10 pe, +0.11 se, +0.18 pb, +0.14 ont).
That gain is the room the strain alleles have, which use the same genome split.

## 3. How protal predicts the unknown share

Source: [Composition.h](../../../src/Profiling/Composition.h).

1. **The bases read.** The aligner counts every read's bases into the SAM header (`@CO protal scanned reads`). A pair's
   overlapping mates count once: the bases are scaled by the fragment base share.
2. **What the called species explain.** Each called species contributes its depth times its genome size:
   - its depth is the fragment bases on its marker genes per marker base (`Taxon::VerticalCoverage`);
   - its genome size is the CheckM-corrected mean of its GTDB genomes (`species_priors.tsv`).

   Their sum over the scanned bases is the explained share.
3. **The unknown genomes.** The bases left unexplained (never below 0) are divided by the called species' average
   genome size. That average is weighted by depth, so it is per cell.
4. **`?`** is the unknown genome equivalents over all genome equivalents. The called species count as the sum of their
   depths. So `?` is a share of cells, not of reads, and it ends the profile; the called species' abundances are
   scaled to the rest.
5. **The species missing** are the unknown genome equivalents over the called species' median depth (and over the
   lowest depth).

## 4. What moves it

`?` is a remainder, so every error in the called species' bases lands in it with the opposite sign.

| Factor | Direction | Size in v19 (median, test and training) |
|---|---|---|
| Depth of the called species, deep samples | pe/se −3% to −8% (100 bp reads lose most) → `?` +0.01-0.02 | pe −3.5%, se −0.5% (console); per point −3% to −10% |
| Shallow samples: called species are the lucky ones | depth overestimated → explained too high, `?` too low | pe 1,000 pairs: depth +100% to +140%, `?` −0.27; 2,000: +50% to +60%, `?` −0.11 |
| PacBio depth (open, below) | +5% at 6 Gb, +12-35% at 25-250 Mb → `?` −0.11 to −0.15 at every depth | pb `?` −0.16 test |
| Nanopore depth | −11% to −19% from 150 Mb up → `?` +0.04 to +0.10 | ont deep samples |
| Calls (the knob) | a missed species' reads go to `?` (by design); an FP species takes reads from `?` | truth "given the calls" separates this out |
| Genome size from the species' GTDB mean | a strain's genome differs: +2% median (10-90%: −3% to +19%) | moves explained by the same % |
| Average genome size used for the unknown | the unknown organisms are assumed to be the called species' size; AGS +5% here | `?` in cells, off where the missing are larger or smaller (soil) |
| Host, eukaryotes, viruses, plasmids | count as unknown bacterial genomes | host samples: `?` +0.85 to +0.92 |
| Reads of missing relatives on a called congener | inflate the congener's depth | no trend with the share the database lacks (table below) |

`?` error (estimate − truth given the calls), samples without a host, by fragments scanned; the truth is 0.36-0.62
here because the training database lacks about 40% of each sample's species:

| fragments | pe | se | pb | ont |
|---|---|---|---|---|
| < 1,500 | −0.264 | −0.273 | −0.316 | −0.204 |
| 1,500-5,000 | −0.112 | −0.110 | −0.241 | −0.030 |
| 5,000-15,000 | −0.040 | −0.045 | −0.180 | |
| 15,000-50,000 | −0.042 | −0.036 | −0.146 | +0.042 |
| 50,000-500,000 | +0.011 | +0.023 | −0.130 | +0.039 |
| > 500,000 | +0.004 | +0.012 | −0.113 (5) | +0.096 (16) |

The species' depth error by the share of the sample's species the training database lacks, samples of 10,000
fragments or more:

| share lacking | pe | pb | ont |
|---|---|---|---|
| 0-0.2 | −0.054 | +0.133 | −0.159 |
| 0.2-0.4 | −0.026 | +0.132 | −0.138 |
| 0.4-0.6 | −0.042 | +0.251 | −0.060 |
| 0.6-0.8 | −0.053 | +0.166 | −0.124 |

**PacBio's bias is not settled.** It shrinks with depth: +84× at 150 kb of reads, +3.4× at 1 Mb, +0.35 at 25 Mb,
+0.12 at 150 Mb, +0.05-0.08 at 1.5-6 Gb. So part of it is the shallow-sample selection, stronger for long reads
because a species has fewer, larger reads on its marker genes. But the deep scenario samples keep +0.13-0.33 in depth
and +0.15-0.27 in explained share (gut, soil, moderate at 1.5-6 Gb), and Nanopore, from the same simulator, errs the
other way.

Two things look ruled out:
- the simulator's placement at contig ends, since it redraws a read that runs off its contig (`LongReadSimulator.cpp`
  `kPlacements`) and would bias both read types alike;
- the reads of the missing species (table above).

Finding the cause needs one PacBio sample's per-species depths against its manifest. The build's profiles stayed on
node-local scratch.

**The missing-species count is rough.** It is 2-3 times the species present and not called (median ratio +1.9). The
unknown share is dominated by abundant species the database lacks, and dividing their genome equivalents by the median
depth counts one abundant species as several.

## 5. Next

- **Pull `504de80` onto the HPC checkout before the next build.** v19 tested nothing new for the classifier.
- **Gut:** read v19's 50 gut FN from `model_logs/error_reads` before blaming the clouds' hold-out.
- **Ancestry:** fixed-site features from `strain_alleles.tsv`, and realistic in-silico strains (section 7).
- **Composition:**
  - a shrinkage of a few-fragment species' depth, or `?` defined as the organisms the database lacks (section 6);
  - one PacBio and one Nanopore sample profiled locally against their manifest, per species, to find the depth biases;
  - then a per-read-type depth calibration from `composition_accuracy*.tsv` if no cause is found;
  - a shallow-sample warning (fewer than ~5,000 fragments for short reads) on `?`;
  - the host's reads left out of `?` when a host is screened.

## 6. Follow-up: shallow samples, counts against calls (2026-10-09)

**The question (the user's).** At shallow depth, does protal take a few-fragment species' abundance from a median
rather than its counts? Would the raw fragment counts be a less biased basis for `?`?

**protal already counts there.** A taxon's depth (`Taxon::VerticalCoverage`, `BlendedDepth` in `Profiler.h`) blends
two estimators:
- its own reads' aligned bases over the length of all its expected marker genes, the genes without reads included;
- the median over the genes hit.

The median gets weight only from 80% of the genes hit (full at 95%), or from a median depth of 0.5-1× with at least a
quarter of the genes hit. A species with a few fragments on 2 of 120 genes gets the count-based depth alone, which is
unbiased for it.

**The total is right; its attribution is not.** The same samples are compared with two truths:
- **the database's split:** the share of the reads from the species the database has, and of the cells from the
  species it lacks;
- **given the calls:** the same, with only the species called and present counted as explained.

| pe, no host, fragments | species (lacking) | called | explained − DB species' share | explained − given the calls | `?` − cells the DB lacks | `?` − given the calls |
|---|---|---|---|---|---|---|
| < 1,500 | 129 (42) | 14 | +0.019 | +0.279 | −0.034 | −0.264 |
| 1,500-5,000 | 115 (39) | 21 | −0.013 | +0.131 | +0.008 | −0.112 |
| 5,000-15,000 | 136 (49) | 37 | −0.026 | +0.058 | +0.028 | −0.040 |
| 15,000-50,000 | 100 (35) | 47 | −0.011 | +0.011 | −0.024 | −0.042 |
| 50,000-500,000 | 114 (39) | 63 | −0.027 | −0.012 | +0.023 | +0.011 |
| > 500,000 | 1667 (874) | 685 | −0.022 | −0.014 | +0.010 | +0.004 |

- **Below 1,500 fragments, 14 called species carry the reads of about 87 database species present.** Their depths are
  extrapolated from the marker genes to the whole genome. So they also account for the reads of the ~73 species with
  too few marker fragments to be called.
- **The sum is right.** The explained share is within 0.02-0.03 of the database's species at every depth. `?`
  therefore estimates, within 0.03 at any depth, the cells of the species the database lacks. It does not estimate
  what the profile leaves unnamed.
- **The error is per species.** At shallow depth, a species is called because its marker count drew high. Its
  count-based depth is the right estimate for a random species, but too high for a called one. That is selection, not
  the estimator.

**So counts alone do not lower the bias; a shrinkage would.** Two ways to get the called species' share right:
- **Shrink a few-fragment species' depth** towards the sample's abundance distribution (empirical Bayes, with Poisson
  counts on its marker length).
  - Fit the prior from the sample's well-covered species, or from all taxa with reads.
  - The posterior mean replaces n/m for small n.
  - The depth taken away is the expected share of the database's species below detection. It can go into `?`, or into
    a line of its own ("known species below detection" beside "unknown").
- **Or keep the depths** and define `?` as the cells of the organisms the database lacks, which is what it measures.
  Then say that the called abundances of a shallow profile include the species below detection.

PacBio is not this case: it is off against the database's species too (+0.12 to +0.26 at every depth), so its depths
are too high.

The raw counts per species are in the outputs:
- `.profile.log` gives every taxon's fragments;
- the training tables give `fragments`, `fragments_all`, `em_fragments` and `depth` per taxon.

The truth per species (the manifests) stayed on the build's scratch. A per-species check of a shrinkage needs a local
collection of shallow samples.

## 7. Follow-up: the ancestry sites of the real strains alone (2026-10-09)

[fixed_sites_check.py](fixed_sites_check.py) ([output](fixed_sites_check.txt)) is the v18 report's script, run on v19.
- It splits the FN taxa's own records by what was simulated (the tables' `meta_rep_genome`, `meta_insilico_strain`).
- It compares each subset with the FP taxa of the same allele status.
- The build's console AUC (0.58-0.68) mixes all three FN classes.
- The alleles are real genomes in both the report and protal: the full reference holds GTDB's genomes only, and an
  in-silico strain's species has one genome.

AUC per taxon (≥ 10 sites) against the FP genus taxa, FN high:

| FN class (species with alleles) | taxa pe | congener sites pe / se / pb / ont | fixed sites pe / se / pb / ont |
|---|---|---|---|
| real strain (yes) | 1,922 | 0.450 / 0.463 / 0.373 / 0.376 | **0.647 / 0.667 / 0.614 / 0.585** |
| real strain (no) | 837 | 0.520 / 0.537 / 0.351 / 0.409 | the same (no polymorphic sites) |
| in-silico strain (no) | 525 | 0.728 / 0.765 / 0.706 / 0.755 | the same |
| representative (any) | 471 | 0.889 / 0.909 / 0.943 / 0.870 | 0.853 / 0.881 / 0.912 / 0.828 |

- **For real strains with alleles, the fixed sites turn the sign.**
  - At all congener sites, the missed strains carry the congener's base more often than the novel congeners do (AUC
    0.37-0.46).
  - At the sites fixed within the species, they carry it less (0.59-0.67).
  - 40-47% of their mismatches fall on the 2.5% of bases that are polymorphic (26% for ont).
  - Without their own genome, about half of the lift v18 reported survives in the hardest pair, FN against FP (+0.20 against
    +0.37 in pe).
- **The in-silico strains are too easy here.**
  - At congener sites they lack the species' base at 3.5% of sites (pe; table in the output). The real strains
    lack it at 10.0%, the FP at 8.7% and the representatives at 0.5%.
  - Their differences are drawn by gene rate and codon (`insilico_strains.py`: conservation factors, kappa,
    omega), not by where the genus varies, so they seldom land where the congeners differ. Real strains share
    the congeners' base there: ancestral polymorphism, recombination, fast sites.
  - Their AUC (0.71-0.77) means the model learns from them that a present strain rarely sides with a congener. Real
    strains do.
  - They are 18% of the simulated non-representative genomes and 10% of the FN records.

**What reaches the classifier.**
- **v19:** none of it. protal's `ancestry` features (`AncestrySites.h`) count every congener site; the polymorphic
  mask exists only in the report.
- **`504de80`:** the species' alleles reach the classifier as `allele_explained_share` and `allele_identity_gain` (the
  share of a read's differences a known strain explains, and the identity it gains), and through the allele-aware
  candidate scores. That covers the polymorphic sites' half of the signal. It does not ask the fixed-site question:
  where the species does not vary, does the read side with the congener?

**What would take full advantage**, in order:
1. **Fixed-site ancestry features in protal**, from `strain_alleles.tsv` (no new table).
   - A site is polymorphic where any stored allele of the copy has an edit there (the edits are sorted by position;
     at most 4 alleles per copy).
   - Count the fixed sites apart, and add `ancestry_fixed_congener_share` (the share of fixed sites with the
     congener's base) and its sites per record.
   - A species without alleles gets 0, as in the alleles group, so that the feature does not say which species have
     other genomes (the cluster-size prior).
   - Better still, express it as a gain: the congener share at all sites minus that at the fixed sites, which is 0
     without alleles.
2. **A three-way site.** Where a stored allele carries the congener's base, a read with that base is explained by a
   known strain, not by the congener. Count it as agreement (a shared polymorphism).
3. **Weights:**
   - weigh a fixed site by how many alleles cover it: 0 of 4 says little, 4 of 4 a lot;
   - sum the sites over a taxon's records instead of averaging per record, so that long reads weigh by their sites;
   - in training, weigh the in-silico strains' rows down (e.g. 0.5) until their sites are realistic.
4. **Realistic in-silico strains** (`insilico_strains.py`).
   - Draw a share of their substitutions at the sites where the genus's copies differ, sometimes to a congener's base.
   - Calibrate on the real strains: the species' base missing at 10% of congener sites, and 40-47% of all their
     mismatches on polymorphic sites.

**Expected size.** protal's mask will know fewer polymorphic sites than the report's, so its gain will be below
0.59-0.67:
- protal's alleles come from half of each species' genomes (the hash split), at most 4 per copy;
- the report's came from up to 6 of all of them, with the own genome left out.

A real strain of a species with one genome gains nothing in either case (837 of 2,759 FN taxa, pe).
