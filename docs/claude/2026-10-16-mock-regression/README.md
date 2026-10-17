# Why protal 0.8.5 scores below 0.6.0a on some real mocks: thin calls of foreign DNA, and divergent strains

The user's mock-community runs on the cluster (`score_mocks.py score --store_results_tar`) show 0.8.5 behind 0.6.0a.
The in-silico test samples say the opposite. The loss has two causes, and neither occurs in the simulations.

- **Thin calls of foreign DNA (meslier: −0.22 F1 on pe).**
  - **What 0.8.5 calls:** species from a handful of fragments. On the Meslier mocks' paired-end samples it calls a
    median of 52 false species per sample, where 0.6.0a calls 4. Recall is the same (0.943).
  - **The false species are real foreign DNA, not misplaced reads.** 146 of the 205 are of genera the mock does not
    have (human-gut Alistipes, Agathobacter, Dialister, Gemmiger; Thermus, Limnothrix). 192 of 205 have reads at least
    97% identical to their reference. Their median is 6 fragments, against 5,718 for the true members. No true member
    of these samples has fewer than 11 fragments.
  - **It depends on the library:** DNBSEQ samples get 52-73 false species, Ion Torrent single-end 116-150, and
    Illumina HiSeq 9-12, from the same mock DNA. The DNBSEQ libraries have a third to a half more taxa with reads than
    the HiSeq ones.
  - **0.6.0a sees the same DNA:** nearly every false species of 0.6.0a is also one of 0.8.5's (28 of 29 in the four
    paired-end samples compared). 0.6.0a stops at the stronger ones; 0.8.5 places about 1.8 times as many fragments on a
    species, and its model accepts 2-10 perfect reads.
  - **The second profiler reports 2-3 false species per sample:** two contaminants of the mock DNA itself
    (*Limnothrix sp001693275*, *Stenotrophomonas maltophilia_S*) that every tool calls, and none of the thin ones.
- **Why the simulations cannot show it:** every read of a simulated sample comes from a community member, one of its
  in-silico strains, a held-out relative or the host. A taxon with a few high-identity reads is therefore often truly
  present there. In the v18 training rows, 32-40% of such rows with 3-10 fragments are present, and 22% of those with
  1-2. Real libraries carry low-level foreign DNA (other samples on the flow cell, reagents, the lab). The model has
  never seen a negative example of it.
- **Divergent strains (treichel: −0.024 F1).**
  - **What 0.8.5 loses:** true species that 0.6.0a called, 117 times, against 65 gains. Mostly *Blautia_A
    fusiformis* (lost in 33 samples), *Bacteroides ovatus* (31) and *Bifidobacterium pseudocatenulatum* (21).
  - **These are not thin:** a median of 151 fragments, reads 98% identical to the reference (median 0.980), MAPQ 21
    (a congener takes a share of each read), and probabilities 0.04-0.69.
  - **Why the simulations miss it:** this is the strain weakness the simulations understate. The Treichel isolates
    are further from GTDB's genomes than the in-silico strains (`--omega`,
    [2026-10-07-error-read-signatures](../2026-10-07-error-read-signatures/README.md)).
- **Elsewhere 0.8.5 is better:** D6331 (0.900 against 0.877) and hCom2 (0.694 against 0.663, with 99 true members called
  on 1-2 fragments).
- **What a floor would do:** dropping 0.8.5's calls below 3·10⁻⁵ of the profile lifts meslier pe from 0.717 to 0.957
  (above 0.6.0a's 0.941), meslier se from 0.513 to 0.898, and D6331 ONT from 0.778 to 0.947. hCom2 falls from 0.694 to
  0.671, and treichel is unchanged (section 5). A fixed minimum of fragments does as well on meslier but costs hCom2 and
  treichel more.

Nothing in protal was changed for this report. Section 6 lists what would settle the open questions.

## Data and commands

- **The archive:** `local/mock_results_2026-10-16.tar.gz` (101 MB), packed on the cluster by `score_mocks.py`
  (commit da69b40) and unpacked to `local/mock_results_2026-10-16/`.
- **The runs** (`README.txt` of the archive):

  | Run | Version | Database | Samples |
  |---|---|---|---|
  | `0.6.0a` | 0.6.0a (paired-end only) | `protal_0.5.1_r226` | 367 |
  | `0.8.5_956e18b` | 0.8.5, commit 956e18b | r226 v31 (gbm models, gates pe 4, se 4, ONT 8, PacBio 16) | 429 |
  | `0.8.4_55f850d` | 0.8.4 | r226 v29 | 7 |
  | `sylph_1.0.0_gtdb-r226-c1000-dbv1` | the second profiler | GTDB r226 | 429 |

  `0.8.5_956e18b_v0` was scored, but none of its files were packed. This is most likely because its set folders are
  links: the packing did not follow links, and scoring did. That is fixed in `score_mocks.py` with this report. Its
  rows are left out here.
- **None of these runs gave protal the truth,** so they have no `.truth_annotated` with the model's features. The
  evidence below comes from `.profile.log` (the Summary column: genes, fragments as "Total Hits", mean ANI, mean MAPQ,
  depth). 0.6.0a's `.profile.log` has no header, so the archive kept only those up to 1 MB (most hCom2 and some
  treichel samples lack them).
- **Scripts** (WSL, Python 3.12, `taskset -c 0-3`):
  - [evidence.py](evidence.py) `local/mock_results_2026-10-16 local/mock_regression`: `calls.tsv`, each run's species
    of every sample with its status and evidence, and for a false species whether its genus is in the sample's truth
    and whether it is a member of another truth table;
  - [bins.py](bins.py) `local/mock_regression/calls.tsv local/mock_regression`: true and false calls by fragments
    (`bins.tsv`), and the median F1 with a minimum of fragments (`min_fragments.tsv`) or of profile share
    (`min_share.tsv`);
  - [thin_rows.py](thin_rows.py) `local/v18/{training,test}/training_data.tsv`: the simulated rows by fragments and top
    identity, with the share present and called. The v31 tables were removed after its report, so these are the r226
    v18 build's paired-end rows (2026-10-08).

## 1. Where 0.8.5 loses

Median per sample (`scores/summary_by_run.tsv`), on the samples each run profiled. 0.6.0a profiled the paired-end
samples only.

| Set | Read type | Samples | F1 0.6.0a | F1 0.8.5 | Second profiler | FP 0.6.0a | FP 0.8.5 | FN 0.6.0a | FN 0.8.5 |
|---|---|---|---|---|---|---|---|---|---|
| d6331_dna | pe | 2 | 0.877 | **0.900** | 0.889 | 1 | 1.5 | 2.5 | 1.5 |
| meslier | pe | 5 | **0.941** | 0.717 | 0.964 | 4 | 52 | 4 | 4 |
| treichel | pe | 84 | **0.933** | 0.909 | 0.933 | 0 | 0.5 | 3 | 4 |
| hcom2_invitro | pe | 276 | 0.663 | **0.694** | 0.655 | 0 | 1 | 58.5 | 54 |
| meslier | se | 5 | | 0.513 | 0.930 | | 124 | | 3 |
| meslier | pb | 1 | | 0.778 | 0.834 | | 32 | | 5 |
| meslier | ont | 2 | | 0.952 | 0.928 | | 3.5 | | 3 |
| d6331_ont | ont | 2 | | 0.778 | 0.871 | | 7 | | 1 |

## 2. The meslier false species follow the sequencing platform

Per sample (`scores/samples_by_run.tsv`), 20M pairs or reads each for the short reads:

| Sample | 0.6.0a FP | 0.8.5 FP | Second profiler FP | Taxa with reads, 0.6.0a | Taxa with reads, 0.8.5 |
|---|---|---|---|---|---|
| mock1 DNBSEQ G400 | 4 | 73 | 2 | 423 | 844 |
| mock1 DNBSEQ T7 | 9 | 52 | 3 | 470 | 900 |
| mock2 DNBSEQ T7 | 13 | 59 | 3 | | 992 |
| mock1 Illumina HiSeq 3000 | 3 | 9 | 2 | 318 | 631 |
| mock2 Illumina HiSeq 3000 | 4 | 12 | 3 | 392 | 718 |
| mock1 Ion Proton / S5 | | 116 / 124 | 6 / 7 | | 838 / 804 |
| mock2 Ion Proton | | 57 | 3 | | 748 |
| mock3 Ion Proton / S5 | | 150 / 132 | 7 / 8 | | 762 / 811 |
| mock1 PacBio Sequel II | | 32 | 18 | | 204 |
| mock1 / mock3 MinION R9 | | 4 / 3 | 0 / 1 | | 187 / 170 |

- **The same mock DNA gives more foreign taxa on DNBSEQ and Ion Torrent:** a third to a half more taxa with reads than
  HiSeq, in either protal version.
- **The false species of 0.6.0a are the stronger part of 0.8.5's:** 8 of 9, 4 of 4, 13 of 13 and 3 of 3 in the four
  paired-end samples compared.
- **The second profiler's false species:** *Limnothrix sp001693275* and *Stenotrophomonas maltophilia_S* in nearly
  every sample, on every platform, called by both protal versions too: contaminants of the mock DNA. Beyond those,
  a few of protal's stronger ones on Ion (all 7 of its false species of mock1 Ion S5 are among 0.8.5's 124).

## 3. The false species are foreign DNA with few reads

True and false calls of the meslier samples by their fragments (`bins.tsv`; summed over the samples):

| Run, reads | Status | 1-2 | 3-5 | 6-10 | 11-25 | 26-100 | 101-1000 | >1000 |
|---|---|---|---|---|---|---|---|---|
| 0.6.0a, pe | TP | 0 | 0 | 0 | 7 | 24 | 47 | 201 |
| 0.6.0a, pe | FP | 0 | 4 | 7 | 0 | 5 | 4 | 0 |
| 0.8.5, pe | TP | 0 | 0 | 0 | 3 | 27 | 72 | 260 |
| 0.8.5, pe | FP | **36** | **61** | **55** | 40 | 5 | 8 | 0 |
| 0.8.5, se | TP | 0 | 1 | 2 | 10 | 26 | 65 | 232 |
| 0.8.5, se | FP | **154** | **177** | **111** | 82 | 48 | 5 | 2 |

Where they come from (0.8.5, meslier):

| Reads | False species | Genus in the mock (a relative's reads) | Member of another truth table | Foreign | Reads ≥ 97% identical |
|---|---|---|---|---|---|
| pe | 205 | 16 | 43 | 146 | 192 |
| se | 579 | 37 | 54 | 488 | 578 |

The 43 members of another truth table are 29 distinct species. Of these, 2 are members of another Meslier mock; the
rest are members of the gut mocks of other sets (Treichel, hCom2, NIBSC), so gut DNA as well.

An example, `meslier_mock1_DNBSEQ_T7` in 0.8.5 (probability, genes, fragments, mean identity):

| Species | P | Genes | Fragments | Identity |
|---|---|---|---|---|
| *Methanothermobacter thermautotrophicus* | 0.72 | 1 | 2 | 1.000 |
| *Cytobacillus pseudoceanisediminis* | 0.77 | 1 | 2 | 1.000 |
| *Thermus thermophilus* | 0.97 | 2 | 4 | 0.998 |
| *Dialister invisus* | 0.96 | 5 | 9 | 0.977 |
| *Agathobacter rectalis* | 0.99 | 9 | 15 | 0.980 |
| *Alistipes putredinis* | 0.81 | 9 | 18 | 0.976 |

Mock 1 has no member of these genera. A read 97-100% identical to *Alistipes* cannot come from its archaea and
thermophiles: these are reads of organisms the sample really contains, at 10⁻⁶ to 10⁻⁵ of its cells.

0.8.5 also finds more reads everywhere. On the true species both versions call in the same sample, 0.8.5 has a median
1.8 times the fragments of 0.6.0a in meslier (10-90%: 1.65-2.10), 1.8 in hCom2, and 1.36 in treichel. It also lists
twice as many taxa with reads (section 2). That is the alignment's gain in sensitivity, and it gives the model about
twice the evidence on the same foreign DNA. (D6331 DNA gave 24 times, which this archive cannot explain; possibly its
0.6.0a run read fewer of the sample's reads.)

## 4. Why the simulations reward thin calls

The training and test rows of a build (here r226 v18, paired-end) by fragments and top identity:

| Fragments | Rows (identity ≥ 0.97) | Present | Called | Rows (< 0.97) | Present |
|---|---|---|---|---|---|
| 1-2 | 55,427 | 0.220 | 0.003 | 176,785 | 0.003 |
| 3-5 | 29,861 | 0.316 | 0.229 | 26,311 | 0.001 |
| 6-10 | 22,826 | 0.397 | 0.395 | 8,588 | 0.002 |
| 11-25 | 26,302 | 0.497 | 0.563 | 4,174 | 0.001 |
| 26-100 | 26,509 | 0.636 | 0.749 | 1,765 | 0.002 |

(The test rows and the deep samples alone, ≥ 5M pairs, give the same picture: 16-20% present at 1-2 fragments, 26-29%
at 3-5.)

- **A thin row is present a fifth to two fifths of the time in simulation.** A row of 3-10 fragments at 97% identity
  or more is present in 32-40% of cases, and the model calls 23-40% of them.
- **The reads of an absent taxon in simulation come from a relative:** their identity is below 0.97 (0.1-0.3%
  present). High-identity reads of an absent taxon barely occur, because no simulated sample holds foreign DNA.
- **The meslier false species are exactly the rows this teaches the model to call:** 3-25 fragments at 97-100%
  identity. In a real library they are almost never members.

These are v18's rows. The v31 build behind 0.8.5 trains on scenarios and pushed hard for thin true species
([2026-10-12 v24b errors](../2026-10-12-r226-v24b-errors/README.md): 60-76% of the misses had fewer than 3 own
reads). That
can only have raised the called share at 1-10 fragments; its tables were removed after its report, so this is not
measured. The knob is no factor: one depth point per model, so a constant 0.69 for pe and 0.59 for se.

## 5. Treichel: divergent strains rejected

Of the species and samples where both versions have a status, 3,233 are TP in both and 309 FN in both. 117 go from TP
in 0.6.0a to FN in 0.8.5, and 65 the other way. The 117 lost have a median of 151 fragments in 0.8.5, identity 0.980
and MAPQ 21. The largest losses (0.8.5's evidence in the samples where it lost them):

| Species | Samples lost | Fragments | Identity | MAPQ | P |
|---|---|---|---|---|---|
| *Blautia_A fusiformis* | 33 | 17-1,362 | 0.969-0.984 | 4-22 | 0.04-0.59 |
| *Bacteroides ovatus* | 31 | 31-528 | 0.974-0.991 | 7-29 | 0.04-0.66 |
| *Bifidobacterium pseudocatenulatum* | 21 | 22-660 | 0.972-0.990 | 20-51 | 0.08-0.68 |
| *Phocaeicola vulgatus* | 9 | 59-253 | 0.980-0.994 | 13-26 | 0.06-0.69 |
| *Oliverpabstia intestinalis* | 7 | 89-1,514 | 0.976-0.983 | 22-40 | 0.06-0.64 |

- **0.6.0a called the same species on fewer fragments at the same identity:** e.g. *B. ovatus* from 30-700
  fragments, identity 0.976-0.995, P 0.54-0.99.
- **The isolates are not GTDB's genomes:** their reads sit about 2% from the reference, and a congener takes a share
  (MAPQ around 20).
- **The 0.8.5 model reads that pattern as a relative's reads:** in simulation it usually is one, an absent congener
  of a present species. The in-silico strains, which should teach the model this case, were found too close to their
  references before (`--omega`, [2026-10-07-error-read-signatures](../2026-10-07-error-read-signatures/README.md)).
  Which features reject them needs the feature tables (section 6).

A share floor leaves treichel as it is (0.909 at every floor in `min_share.tsv`): its losses are misses, not thin
false calls.

## 6. What a floor would do, and what would settle it

Median F1 of 0.8.5 if the calls below a share of the profile were dropped (`min_share.tsv`; a dropped true call
counts as missed; bold: 0.8.5's best floor for the set):

| Set, reads | No floor | 10⁻⁵ | 3·10⁻⁵ | 10⁻⁴ | 0.6.0a | Second profiler |
|---|---|---|---|---|---|---|
| meslier pe | 0.717 | 0.936 | **0.957** | 0.940 | 0.941 | 0.964 |
| meslier se | 0.513 | 0.777 | 0.898 | **0.936** | | 0.930 |
| meslier pb | 0.778 | 0.815 | **0.851** | 0.837 | | 0.834 |
| meslier ont | **0.952** | 0.951 | 0.946 | 0.945 | | 0.928 |
| d6331_ont | 0.778 | 0.824 | **0.947** | 0.929 | | 0.871 |
| d6331_dna pe | 0.900 | 0.900 | **0.915** | 0.897 | 0.877 | 0.889 |
| hcom2 pe | **0.694** | 0.682 | 0.671 | 0.640 | 0.663 | 0.655 |
| treichel pe | 0.909 | 0.909 | 0.909 | 0.909 | 0.933 | 0.933 |

A minimum of fragments (`min_fragments.tsv`) does as well on meslier pe (26 fragments: 0.957). But it costs hCom2 more
(0.663) and treichel too (0.893). The fragments a species needs grow with the sample's depth, which a share follows.

What would settle the open questions:

1. **The features:** profile the meslier samples (17) and treichel's paired-end samples again with the current
   `score_mocks.py run`. It gives protal the truth, so each sample gets `.truth_annotated` with every taxon's features.
   Then [swap_features.py](../2026-10-15-mock-near-misses/swap_features.py) can tell which features carry the thin
   foreign calls and the rejected strains, as for the D6331 near-misses.
2. **Foreign DNA in training:** simulate it. Add to each simulated sample a few reads (Poisson, a few to a few dozen
   fragments a taxon) of a few unrelated genomes: gut species, reagent genera such as *Thermus* and *Ralstonia*, and
   other samples' species, labelled absent. This is the principled fix: the model then learns what a few perfect reads
   of an unrelated genus mean in a real library. The test set should carry it too, or simulated F1 keeps rewarding
   thin calls.
   *2026-10-17: decided otherwise.* The user's call is that contaminants are real species in the sample, so they are
   true positives, not absent ones. They were added to the simulations labelled present
   ([2026-10-17-contaminants](../2026-10-17-contaminants/README.md)). On that view, the meslier "false" species of
   section 3 are, for the most part, real organisms the truth tables do not list.
3. **A detection floor as an option** (e.g. `--min_share 3e-5`), until the model knows foreign DNA. It cost hCom2
   0.023 here and gained meslier 0.24-0.39. hCom2's thin true members are real (99 called on 1-2 fragments), so the
   default needs care.
4. **Strains:** real strain genomes, or a larger `--omega`, in the training simulations, for the treichel misses.

## Caveats

- **Few samples:** meslier has 5 paired-end and 5 single-end samples; its medians move with one sample.
- **The truth has no contamination:** the truth tables hold the designed members. A "false" species here is often a
  real organism in the library, and below any detection limit a user would want. The benchmark counts it all the same,
  as it should for a profiler of communities.
- **0.6.0a used an older database** (`protal_0.5.1_r226`), so part of its lower sensitivity is its database, not its
  model. This report does not separate the two.
- **0.6.0a's evidence is partial:** its `.profile.log` files over 1 MB are not in the archive (no header to cut by), so
  its rows in the evidence tables are incomplete for hCom2 and treichel. Its scores are complete.
