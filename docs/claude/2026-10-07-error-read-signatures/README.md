# The reads behind the false positives and false negatives: what their alignments tell

**Question (2026-10-07).** The r226 v15 build kept the SAM records of the reads behind each model's false positives
(FP) and false negatives (FN). Do these reads and their alignments, as protal reports them, carry signatures that
would let the presence model classify FP and FN taxa better?

**Short answer.** Hardly. A read counted for an FP taxon is almost always a read of a species held out of the
training database, landing on its closest congener at about 97% identity. The reads of a missed taxon are mostly
reads of a divergent strain of a wide species, at about 96-97%. As alignments to the marker genes, the two cannot be
told apart, read by read or taxon by taxon:
- tags (ZU, ZT, ZA, ZF), MAPQ, clips, gene ends, mates, gene, position and composition give AUC 0.57-0.64;
- the missed strains are the *more* divergent of the two;
- the features that point the right way describe the species, not the reads: how many genomes its GTDB cluster has,
  how complete its representative is.

The way forward is the database (strain alleles of wide species; a per-species marker radius) and the reporting
(novel congeners as genus-level calls), not more read features.

**Follow-ups (same day).**
- **Sensitivity (section 5).** FP taxa are not a sensitivity problem: 98-99% of their reads come from species the
  database lacks. For FN taxa, alignment loses 2-3% of the reads tried against the own species (Nanopore 15%). More
  is lost before that: as many own reads go to a congener as stay, the own species never tried, and more so the
  larger the genus. That points at the top-3 candidate limit (`--align_top`).
- **Genus-level calls (section 6).** If every current FP were reported at genus level, 73-86% of those calls would be
  right (species F1 +0.009 to +0.022, a ceiling no rule reaches). A learned second model relabels 15-440 calls at
  78-95% precision: species F1 +0.002 to +0.003, mostly in soil. The larger value is *adding* genus calls for genera
  with reads and no species call: 82-92% right on the test sets, finding 68-87% of the genera with a novel congener.
  Where novel congeners are rare (the gut scenario), it drops to 64%.
- **Gene uniqueness (section 7).** The beyond-genus singletons hit copies that are unique within the database. A
  graded uniqueness over the database's marker copies would see few of them. Whole genomes, the build's simulations,
  or a test of each copy against its own congeners' copies (section 9) could see more.

## Data and setup

- **Build.** r226 v15: SLURM 24027266, built from `ce85bd7`. Model logs and tables are in `local/v15`; the
  [v15 report](../2026-10-07-r226-v15/README.md) covers it.
- **The records.** The build's `model_logs/error_reads` on the cluster held the SAMs that `error_reads.py`
  (`652ed53`) wrote: one `<sample>.sam.zst` per sample, every record of each taken fragment, tagged `xg`/`xs`/`xe`.
  This includes pe: its step failed, but only after writing 418 of the samples' SAMs.
  [`extract_error_records.sh`](extract_error_records.sh) filtered them on the cluster:
  - it kept FP and FN records, and `seeded:`/`source:` records only of FN taxa (no reads of unseen species);
  - it set QUAL to `*` and trimmed the `@SQ` lines;
  - the result is `local/v15/v15_error_records.tar.gz` (375 MB), unpacked in WSL `~/v15/v15_error_records`.

  | read type | SAMs | FP taxa | FN taxa | records | counted fragments (FP / FN taxon / other) |
  |---|---|---|---|---|---|
  | pe | 418 | 2,367 | 3,988 | 3,210,669 | 84,515 / 85,756 / 55,155 |
  | se | 405 | 2,304 | 1,811 | 2,158,090 | 125,193 / 70,166 / 22,825 |
  | pb | 262 | 2,482 | 3,004 | 565,565 | 17,224 / 19,300 / 56,217 |
  | ont | 264 | 3,245 | 3,334 | 636,065 | 23,184 / 22,424 / 46,416 |

- **Samples by set and name.** A scenario's hold-in sample (set training) and its hold-out sample (set test) have
  the same name, e.g. `sc_soil_pe_p20000000_s_3`: 31 pe names, 29 se, 9 pb, 9 ont, holding 53-63% of the error taxa.
  The first version of this report keyed samples by name alone, which merged each such pair: FP and FN roles leaked
  between the two, and 3,394 pe read names collided. Every script now keys a sample as `<set>:<sample>`
  (`fragments.load`, `fragments.model_tables`), and all numbers here are from the re-run. They moved by at most a
  few points; no conclusion changed.

- **What is counted.** A fragment (a pe pair, or a whole long read) counts for the taxon of its best record (no flag
  0x904) at MAPQ 4 or more, as the profiler counts it. Its source (`xs`) gives the relation of that taxon to the
  read's species: own, same genus, family, order, or other. `heldout_species.txt` says whether the source was held
  out of the training database; `internal_taxonomy.dmp` says whether its genome is the representative, another genome,
  or an in-silico strain.
- **Fragment groups:**
  - *wrong (FP)*: counted for an FP taxon;
  - *own (FN)*: counted for an FN taxon from that species;
  - *own (TP)*: counted for another taxon from that species. These are reads of called present taxa, but only those
    that also seeded on an FN taxon (`zf_fn` ≈ 1), a biased sample;
  - *wrong (other)*: misassigned reads that made no FP.
- **Weighting.** `error_reads.py` of `652ed53` kept every FP-taxon record, so a few FP taxa with thousands of fragments
  dominate the fragment counts (pe: median 5, maximum 12,586). The per-taxon views and the models capped at 20
  fragments per taxon ([`followup.py`](followup.py)) are the ones to trust.
- **Scripts.** Run in WSL with `~/soil13/venv`, one pinned core, niced. Paths are relative to the repository;
  `~/v15/err_analysis` is in WSL.

      bash docs/claude/2026-10-07-error-read-signatures/extract_error_records.sh MODEL_LOGS/error_reads OUT.tar.gz   # cluster
      python3 parse_error_records.py --dir ~/v15/v15_error_records --out ~/v15/err_analysis   # one table per read type, 8 min
      python3 fragments.py --records ~/v15/err_analysis --build local/v15                      # fragments, 5 min
      python3 signatures.py --records ~/v15/err_analysis --build local/v15 > signatures.txt
      python3 followup.py --records ~/v15/err_analysis --build local/v15 > followup.txt
      python3 followup2.py --records ~/v15/err_analysis --build local/v15 --error-records ~/v15/v15_error_records > followup2.txt
      python3 checks.py --records ~/v15/err_analysis --build local/v15 > checks.txt
      python3 sensitivity.py --records ~/v15/err_analysis --build local/v15 > sensitivity.txt   # section 5, 1 min
      OMP_NUM_THREADS=1 python3 unknown_genus.py --build local/v15 > unknown_genus.txt         # section 6, 25 min
      python3 uniqueness.py --records ~/v15/err_analysis --build local/v15 > uniqueness.txt     # section 7, 1 min

  Outputs: [signatures.txt](signatures.txt), [followup.txt](followup.txt), [followup2.txt](followup2.txt),
  [checks.txt](checks.txt), [sensitivity.txt](sensitivity.txt), [unknown_genus.txt](unknown_genus.txt),
  [uniqueness.txt](uniqueness.txt).

## Results

### 1. False positives: reads of a novel congener

| | pe | se | pb | ont |
|---|---|---|---|---|
| FP taxa with counted fragments | 2,367 | 2,304 | 2,425 | 3,193 |
| fragments per FP taxon, median (90%) | 5 (45) | 15 (111) | 2 (11) | 2 (11) |
| FP taxa whose main source is a congener (same genus) | 0.803 | 0.921 | 0.813 | 0.813 |
| ... beyond the genus (family, order or further) | 0.197 | 0.079 | 0.187 | 0.187 |
| FP taxa fed only by held-out species | 0.842 | 0.817 | 0.890 | 0.887 |
| FP taxa fed by one source species | 0.818 | 0.751 | 0.932 | 0.927 |
| FP taxa fed by an FN species of the same sample | 42 (1.8%) | 62 (2.7%) | 15 (0.6%) | 20 (0.6%) |

- **Congeners dominate.** By fragments, 97-99.6% of the FP reads come from the FP taxon's genus, and 98-99% of all
  FP reads from species held out of the training database. They align at identity 0.973 (pe), 0.976 (se), 0.978 (pb)
  and 0.956 (ont), with MAPQ around 100: the novel species' nearest relative is a unique best hit. Swaps, where an FN
  species' reads make the FP of the same sample, are rare.
- **Beyond-genus fragments are a second, small class.** In each read type, 236-604 FP taxa have a fragment from
  another family or further. It is one fragment per taxon (se: two), at identity 0.98-0.99 (pb 0.994, ont 0.96-0.98),
  and 62-70% of them come from held-out species. The rest of the read does not disagree: a long read's share of
  aligned bases on the FP taxon is 0.98-0.99, as for own reads. The only difference is that such long reads carry
  fewer other genes: 0.47-0.60 supplementary records against 0.83-0.99. These are the reference-contamination and
  transferred-gene class of the [2026-10-03 FP anatomy](../2026-10-03-false-positive-anatomy/README.md); section 7
  asks whether a gene-uniqueness score would catch them.

### 2. False negatives: divergent strains of wide species

- **Mostly strains.** 89-97% of the counted fragments of FN taxa come from a genome other than the representative:
  84-91% from another real genome of the species, 6-7% from an in-silico strain. By main source, the FN taxa are
  66-79% strains, 11-14% in-silico strains and 8-21% representatives.
- **Reads that seeded on an FN taxon and failed there are other species' reads.** Of the 411,000-633,000 fragments
  per read type that seeded on an FN taxon and did not align to it, 99.3-99.9% are other species' reads, about 60% of
  them from held-out species; the FN species' own such reads are a median 1 per taxon (90%: 3-5). FN taxa simply sit
  among novel congeners. `seeded_not_aligned` (median 113 per pe FN taxon, against 3 counted) and
  `failed_candidate_rate` (FN 0.94, FP 0.95, TP 0.84 in pe) therefore say little about presence. Section 5 follows
  the FN species' own reads.
- **Few counted fragments, but not always.** Most FN taxa have few counted fragments: median 3 (pe), 5 (se) and 2 (pb,
  ont); in the pe design, 45% have one. Yet 593 pe FN taxa (15%), 428 se, 148 pb and 166 ont have 20 or more, mostly
  in soil. Their p is a median 0.46 (pe, knob 0.8) and 0.2 for the other types. In only 8-14% of them did their own
  reads go mostly elsewhere, and almost never to an FP taxon.
- **The recurrent ones are wide species.**
  - *JAQVHD01 sp028696375* is missed in 10 pe samples with 20 or more fragments. It has 7 genomes, minimum cluster
    ANI 95.3%; present, it has identity 0.961 and top identity 1.0.
  - *Fimenecus sp004556705* has 31 genomes and minimum ANI 95.8%. Present, it has 55 fragments at 0.965; absent, 20
    fragments at 0.958 from a novel congener.

### 3. Read level: FP reads look like the missed taxa's own reads

Counted fragments, capped at 20 per taxon, without the reads' composition, scored by a boosted model in sample-grouped
5-fold cross-validation:

| | pe | se | pb | ont |
|---|---|---|---|---|
| AUC, FP reads against own (FN) reads | 0.572 | 0.578 | 0.637 | 0.586 |
| identity: FP / own FN / own TP | 0.974 / 0.974 / 0.985 | 0.974 / 0.976 / 0.988 | 0.978 / 0.975 / 0.996 | 0.950 / 0.949 / 0.966 |

- **Against TP reads.** FP reads separate from the biased TP reads well (0.91-1.00), almost all through `zf_fn`, the
  selection itself. Identity and unique k-mers (`ZU`, `ZT`) separate FP from TP, which the model already does.
- **No signal from the other read features** (AUC 0.45-0.55, most within 0.02 of 0.5):
  - `ZA`: alternatives, their least edits more, ties, and negative values (an alternative fits better than the
    alignment the pair or read consensus kept);
  - `ZF`: failed seeds, and failed seeds on an FN or FP taxon;
  - MAPQ, clipping, gene ends, position in the gene, indels;
  - mates: proper pairs, lost or split mates, template length, two genes;
  - long reads: `ZR`, supplementary records, the read's share of aligned bases on the taxon.
- **No FP-prone genes.** The 10 genes richest in FP fragments hold 20% of them in pe and se, and the same 20-21% of
  the own reads; in pb and ont 17-20% against 14-15%. A gene's enrichment correlates only 0.07-0.28 with its
  congeners' identical share.
- **Composition is not a signal.** Fragment-weighted, FP reads have GC 0.425 against 0.50. Capped per taxon this falls
  to AUC 0.43-0.47: it marks which genera have heavy FP taxa, not the errors.

### 4. Taxon level: FP and FN taxa in the model's existing features

Medians of the FP and FN taxa's rows in the training and test tables, against TP and TN rows of the same samples
(pe; the other read types are alike, [followup2.txt](followup2.txt)):

| feature | TP | FP | FN | TN | AUC FN vs FP |
|---|---|---|---|---|---|
| identity | 0.990 | 0.979 | 0.969 | 0.931 | 0.31 |
| lu_per_kb | 82.7 | 47.9 | 23.3 | 0 | 0.27 |
| excess_scaled_median | -0.0003 | 0.0082 | 0.0162 | 0.055 | 0.69 |
| low_mapq_share | 0.018 | 0.122 | 0.449 | 0.400 | 0.61 |
| cluster_genomes_log10 | 0.30 | 0.30 | 0.48 | 0 | 0.61 |
| rep_completeness | 97.4 | 97.2 | 99.8 | 97.3 | 0.59 |
| cluster_ani_radius | 95 | 95 | 95 | 95 | 0.48 |

- **FP taxa look more present than FN taxa on every read-derived feature:** higher identity, more unique k-mers, less
  excess divergence, fewer low-MAPQ reads. Only the species' features point the right way: an FN species' GTDB
  cluster has more genomes (AUC 0.61-0.64 over the read types), and its representative is more complete.
- **Among rows with 20 or more fragments** (pe: 47,710 TP, 586 FN, 398 FP), FN taxa have identity 0.962 against FP
  0.970 and TP 0.990, and `cluster_genomes_log10` 0.60 against 0.30.
- **`cluster_ani_radius` is 95 in every row of every read type**, so the feature carries no information.
- **Two more checks found nothing:**
  - Within- vs between-species gene divergence: the per-gene divergence of FP and FN taxa (8 or more fragments on 4
    or more genes) fits the factors of `gene_congeners.tsv` equally badly (AUC 0.50-0.53). The two factors correlate
    0.76 over the 168 genes.
  - Spread of divergence over genes: at matched mean divergence, AUC 0.42-0.54 for pe/se and 0.22-0.53 for pb/ont. A
    logistic model of divergence and fragments gains nothing from it (pe 0.788 → 0.787), at most 0.012 for PacBio
    (0.733 → 0.745, with the share of genes without a difference).

### 5. Is it sensitivity? The reads the aligner lost

Follow-up question: are these errors reads that seeded but failed to align? [`sensitivity.py`](sensitivity.py)
([sensitivity.txt](sensitivity.txt)) follows every read of the FP taxa back to its species, and every own read of
the FN species to where it went. The v15 SAMs hold all of them: `error_reads.py` of `652ed53` had no cap, and the
counts equal the taxa tables'. A read's ZF tag lists the taxa it was aligned against and failed on.

**FP taxa: no.** 98-99% of their counted fragments come from species the database lacks (pe 0.993, se 0.991, pb
0.981, ont 0.980), and 90-92% of FP taxa are fed mainly by such species. Nothing in the database was missed there:
the congener is the closest reference, and the alignment to it is right. Of the few FP fragments from a species in
the database (0.7-2%), 85-91% never seeded on that species. Across all four read types, only 0-2 FP taxa got most of
their fragments from reads that seeded on their own species and failed there.

**FN taxa: alignment loses little; candidate selection loses more.** The fate of the FN species' own fragments:

| | pe | se | pb | ont |
|---|---|---|---|---|
| own fragments with a record | 2,210,473 | 1,445,049 | 92,393 | 143,234 |
| best on the taxon | 64,773 | 43,141 | 14,528 | 16,335 |
| tried on the taxon and failed: then elsewhere / nowhere | 359 / 847 | 227 / 784 | 87 / 394 | 217 / 2,754 |
| **alignment sensitivity on the taxon** | **0.982** | **0.977** | **0.968** | **0.846** |
| best on another taxon, own species never tried | 74,547 | 37,468 | 8,721 | 7,876 |
| best on another taxon, own species tried and beaten | 289 | 0 | 912 | 723 |
| seeded on unrelated taxa only, aligned nowhere | 2,069,658 | 1,363,429 | 67,751 | 115,329 |

The last row is reads from outside the marker genes with chance seeds: no loss.

- **Alignment.** A read tried against its own species aligns there 97-98% of the time with short reads and PacBio.
  The reads that fail are at 0.93-0.95 identity where they align elsewhere. Nanopore is the exception: 15% fail
  (0.90 identity elsewhere), and 9.6% of ONT FN taxa lose half or more of their reads this way. That is the ONT error
  profile meeting the alignment thresholds, a sensitivity loss of ONT alone.
- **Seeding.** Reads that seed on nothing leave no record. Own fragments per 1,000 simulated read pairs (pe) do not
  fall with divergence: 24 for strains, 15 for representatives, Spearman −0.25 with the reads' identity. Genome sizes
  differ between the groups, but nothing suggests divergent strains lose reads before seeding.
- **Candidate selection.** As many own reads land on a congener as on the species (pe 75,195 against 64,773):
  - they align there at the same identity (0.977 against 0.976 on the species), 97% of them from strains;
  - the own species was never aligned against: it was beaten in only 289 pe reads;
  - protal aligns a read against its top 3 anchors only (`--align_top`, default 3);
  - the loss rises with the size of the taxon's genus in the database: per FN taxon (5 or more aligned own
    fragments), the median share elsewhere is 0 in one-species genera, 0.09-0.12 with 4-10 species, 0.17-0.29 with
    11-30, 0.33-0.40 with 31-100 and 0.25-0.47 above (Spearman 0.26-0.42);
  - 79-89% of the congeners that take these reads are absent and not called: the reads make no FP, they are only lost
    to the FN;
  - per FN taxon, a median 6% (pe) or 22% (se) of the aligned own reads go elsewhere, and 14-30% of FN taxa lose half
    or more.

  The SAMs cannot tell this apart from a strain whose gene is genuinely nearer a congener's copy, which a strain at
  95% ANI can be. A re-run of some soil test samples with `-c 10` would: own reads on the FN taxa, FN count and F1
  against `-c 3`. The ZN tag, which counts the taxa a read's seeds could not tell apart including untried ones,
  would show the cut directly. It came with `11cdf6b`, after v15, so [`crowding.py`](crowding.py) finds none here;
  it reads ZN from the next build's SAMs.

### 6. Genus-level calls for novel congeners: how precise, and do they raise F1?

Follow-up question: how precise would an "unknown species of genus G" call be with the current model, and would it
raise F1? The [2026-10-03 genus fallback](../2026-10-03-genus-fallback/README.md) tried this on v5: a second
forest on the called taxa relabelled 0-1 calls per run, +0.0002 species F1.

[`unknown_genus.py`](unknown_genus.py) ([unknown_genus.txt](unknown_genus.txt)) works on every row of every v15
sample:
- **Scores.** Training rows with species held out (`p_species`, as `calls.tsv.gz`); test rows (the test set and the
  scenarios' hold-out samples) by the final model. Calls at the model's knob: pe 0.8, the others 0.5.
- **Truth.** A genus-level call is one per sample and genus. It is right when the sample holds a species of that
  genus that the database lacks (`meta_novel_congener`).
- **Rules tried:**
  - *oracle*: every FP reported at genus level;
  - *p band* and *identity*: calls below a cut relabelled;
  - *learned*: a boosted classifier of "absent, in a genus with a novel congener", on the rows' 108 features, the
    model's score, and the genus's other rows in the sample. Out-of-fold by sample on the training rows; fitted on all
    of them for the test rows. It *relabels* calls scoring 0.5 or more, or *adds* a genus call for rows not called;
  - *rule*: a genus with reads but no call and at least k fragments gets a genus call.

Test rows, cut 0.5 ([genus_level.py](genus_level.py) in section 10 is an independent cross-check of the relabelling):

| | pe | se | pb | ont |
|---|---|---|---|---|
| species F1 at the knob | 0.9566 | 0.9622 | 0.9785 | 0.9514 |
| FP calls; in a genus with a novel congener | 835; 0.86 | 134; 0.74 | 39; 0.85 | 1,023; 0.87 |
| *oracle*: species F1; genus calls right | 0.9710; 0.85 | 0.9769; 0.73 | 0.9877; 0.85 | 0.9734; 0.86 |
| *learned relabel*: species F1 | **0.9587** | **0.9644** | **0.9801** | **0.9547** |
| ... calls relabelled; right | 351; 0.95 | 60; 0.78 | 15; 0.80 | 440; 0.92 |
| *learned add*: genus calls; right | 15,492; **0.87** | 1,138; **0.82** | 468; **0.91** | 11,074; **0.92** |
| ... genera with a novel congener found | 13,541 of 15,629 | 934 of 1,370 | 426 of 592 | 10,213 of 12,183 |
| *rule*, 10 or more fragments: genus calls; right | 5,296; 0.87 | 504; 0.74 | 147; 0.97 | 1,991; 0.97 |

- **Today's ceiling.** If every FP were reported at genus level instead, 73-86% of those genus calls would be right
  and species F1 would rise by 0.009-0.022. That needs a perfect FP detector, which the features do not give.
- **Simple relabel rules lose.** Relabelling calls of low p or low identity takes more TPs than FPs at every cut:
  species F1 falls, and the genus calls are only 25-75% right.
- **A learned relabel gains a little:** +0.0016 to +0.0033 species F1 on the test rows, relabelling 15-440 calls of
  which 78-95% are right (training rows, out of fold: +0.0015 to +0.0057). The genus-context features do not drive it (pe 0.9585
  without them). The gain sits in the soil scenarios (pe soil 0.9655 → 0.9677, soil_shallow 0.9396 → 0.9425; design
  0.9645 → 0.9648; gut none). That is why it beats v5's +0.0002: v5 had no scenarios. It is a second-stage filter on
  the same rows, so feeding the presence model what the second model sees might give part of it without a second
  model.
- **The larger value is adding genus calls where no species is called.** At cut 0.5, 82-92% of the added genus calls
  are right, and they find 68-87% of the sample-genera that hold a novel congener and have reads. Species F1 is
  unchanged. Counted as positives, the novel congeners would raise an extended F1 from 0.75 to 0.93 (pe test). That
  is a different metric from protal's, though: the novel-congener genera are more than half as many as the present
  species.
- **Where novel congeners are rare, precision falls.** In the gut scenario the added genus calls are 64% right (pe 70
  calls, ont 62), and in the host scenario none of 1-5 calls is. Half to two thirds of the wrong genus calls fall in
  a genus with a present database species (a missed strain, say), where "an organism of G" is still true. The simple
  rule needs the model: in deep gut samples, reads spill into many genera, and "no call, 3 or more fragments" is 4%
  right (pe) or 26% (ont).
- **The simulation sets the prevalence.** The design samples hold novel congeners by construction, and the soil
  scenarios hold hundreds each (pe training: 77,124 novel congeners in 333 samples). A real sample against GTDB r226 has fewer, the more
  so in well-covered communities like the human gut. Precision there would lie nearer the gut scenario's than the
  table's. A real-data check would show it, e.g. how many genus calls a gut metagenome gets.

### 7. Would a gene-uniqueness score catch the beyond-genus singletons?

Follow-up question: should the binary suspect-copy flag become a graded score of how unique each gene copy is in the
database, used to weigh a read's hit?

**How the flag works now.** `--build` compares each species' copy of each marker gene with every other species'
copy, using 12-mer bottom sketches (`GeneIncongruence.h`). It flags a copy within 0.02 of another genus's copy and
at least 0.02 farther from its own congeners'; a run drops records on flagged copies. At r226 that was 5,589 of
14.5M copies (0.04%). The scan already computes, for every copy, the distance to the nearest copy of another genus
and to the nearest congener's copy. It reports the 1.4M cross-genus pairs within 0.05 in `gene_incongruence.tsv`,
but keeps only the flags.

**What the beyond-genus FP reads show** ([`uniqueness.py`](uniqueness.py), [uniqueness.txt](uniqueness.txt)). A copy
that another genus shares near-identically shows in a read's ZA tag: that species aligns too, within a few edits.

| | pe | se | pb | ont |
|---|---|---|---|---|
| beyond-genus FP fragments | 781 | 573 | 538 | 745 |
| with any ZA alternative | 0.19 | 0.10 | 0.02 | 0.02 |
| another genus within 2 edits: these / own reads | 0.024 / 0.005 | 0.009 / 0.001 | 0.006 / 0.001 | 0.008 / 0.001 |
| unique k-mers (ZU) on the taxon, median: these / own reads | 19 / 8 | 12 / 8 | 104 / 113 | 51 / 53 |
| from held-out species | 0.62 | 0.65 | 0.66 | 0.70 |
| on a copy that took such a read in 2 or more samples | 0.23 | 0.36 | 0.35 | 0.32 |
| test set: on a copy that took one in a training sample | 0.16 | 0.40 | 0.23 | 0.21 |

- **Within the database, these copies are unique.** For 98-99% of the reads, no other genus's copy fits within 2
  edits, and the reads carry as many unique k-mers as own reads, or more. At the scan's looser 0.05, 17-22% of them
  have another genus's copy (section 9). A uniqueness score computed over the database's own copies would rate most
  of these hits as fully informative. Section 9's congruence test against the copy's own congeners is the in-database
  check that could reach further.
- **The near-identical sequence is outside the database's marker copies.** 62-70% of these reads come from species
  the database lacks. For the rest, the source species' own copy is not an alternative either (0-2% within 2 edits).
  The read comes from elsewhere in its genome: a second copy, a transferred gene, or the stretch that a contaminating
  contig in the reference was taken from.
- **So a scan would have to look into whole genomes.** For each marker copy, it would ask whether genomes of other
  genera hold a near-identical stretch. The build's gene-neighbour pass reads every genome at hand, but places each
  species' genes only in that species' genomes. Alternatively, a per-copy rate could be learned from the build's own
  simulations: the share of a copy's reads that come from other genera. On the test set, that would have seen 16-40%
  of these reads, a lower bound since only FP taxa's reads are in these SAMs. The class is 8-20% of FP taxa, so at
  most 2-8% of all FP taxa.
- **A graded score in general.** The distances behind it come free from the scan. As a taxon feature (for example,
  the share of a taxon's fragments on copies with another genus's copy within 0.05), it would let the model weigh
  cross-genus copies instead of dropping only the 0.04%. These data do not predict a gain from it on the FP and FN
  taxa here, though. Their reads are not on shared copies, and the read-level version of the same signal
  (`other_genus_fit_share`, `congener_fit_share` from ZA) is already a feature.

### 8. The mutation spectrum: what differs at equal identity

[`calls.py`](calls.py) joins the build's tables with the scores protal calls on (training rows with species held out,
test rows and the scenarios' hold-out samples by the final model, at the knob) and classes each row by the collector's
meta columns: present representative / strain / in-silico strain; absent *novel congener* (`meta_novel_level`
species and `meta_relative_rank` genus); absent *near novel* (the closest novel species at another rank); absent.
[`spectrum.py`](spectrum.py) ([spectrum.txt](spectrum.txt)) looks at `third_position_share`, the share of a
taxon's mismatches at third codon positions (1/3 without mismatches), and the gene-rate features. Medians by class,
20 or more fragments, pe:

| | absent novel congener | absent, near novel | present strain | present in-silico strain | present rep |
|---|---|---|---|---|---|
| identity | 0.941 | 0.931 | 0.988 | 0.989 | 0.991 |
| third_position_share | 0.76 | 0.69 | 0.45 | 0.41 | 0.34 |
| excess_scaled_median | 0.045 | 0.055 | 0.0007 | 0.0006 | -0.0008 |
| lu_per_kb | 13.0 | 7.6 | 77.4 | 85.0 | 92.2 |
| em_own_share | 0.97 | 0.99 | 1.00 | 1.00 | 1.00 |
| low_mapq_share | 0.33 | 0.27 | 0.05 | 0.02 | 0.03 |

- **At the errors it points the other way from the model's reading of it.** Among FN and FP taxa,
  `third_position_share` is *higher* in the FN (AUC 0.66-0.69 pe, 0.55 se, 0.51-0.53 pb, 0.69-0.71 ont, rows with 5
  and with 20 or more fragments), as `excess_scaled_median` is (0.71-0.82): the missed strains are the more divergent
  organisms, and their mismatches are the real, synonymous-rich kind, while a representative's mismatches are
  sequencing errors on all three positions alike. The model has learned "high third-position share → absent" from
  the population, where it holds; at the boundary it does not.
- **Inside the identity band 0.95-0.985**, where strains and novel congeners overlap, the features that still
  separate present strains from absent novel congeners (AUC): identity 0.85-0.93, `lu_per_kb` 0.80-0.88,
  `em_own_share` 0.63-0.79, `cluster_genomes_log10` 0.70-0.74; `third_position_share` 0.80-0.84 (pe), 0.72 (se), 0.53
  (pb), 0.93-0.94 (ont), the novel congeners high.
- **But it adds nothing for real strains** ([`band.py`](band.py), [band.txt](band.txt)). The band holds 70% of
  pe's FP and 71% of its FN (se 74% and 68%, pb 46% and 62%, ont 48% and 34%). In narrow identity bins (0.01 wide)
  the medians of real strains and novel congeners differ by 0.02-0.05 only (pe 0.72/0.70/0.63/0.57 against
  0.74/0.71/0.67/0.62 from 0.95 up; pb 0.82-0.84 against 0.84; ONT's strains sit at 0.37-0.47 against 0.65-0.70
  because their mismatches are errors), and a boosted classifier of present-against-novel-congener inside the band,
  samples grouped, gains from the feature only through the in-silico strains: with every model feature, AUC pe
  0.9885 → 0.9887, se 0.9902 → 0.9909, pb 0.9707 → 0.9737, ont 0.9930 → 0.9936; with the in-silico strains left out,
  base features 0.9798 → 0.9802 (pe), 0.9832 → 0.9831 (se), 0.9458 → 0.9455 (pb), 0.9885 → 0.9889 (ont). What
  separates the band for real strains is unique k-mers (`lu_per_kb` 0.73-0.90 within bins), the EM share
  (0.62-0.89) and the cluster's size (0.60-0.74); identity within a bin 0.5-0.7. So the spectrum is not a lever, and
  the pb "gain" is the model telling in-silico strains from everything else by their spectrum, which no real organism
  has: the next point.
- **The in-silico strains have another spectrum than the real ones.** At 0.95-0.97 identity (10 or more fragments)
  real strains show a third-position share of 0.70 (pe), 0.78 (se), 0.83 (pb); the in-silico strains 0.62, 0.67,
  0.79; at 0.97-0.98, 0.63/0.76/0.84 against 0.56/0.66/0.72. `insilico_strains.py` keeps an amino-acid change with
  probability `--omega` 0.15 and draws per-frame rates from a gamma; the real strains' differences are more
  synonymous than that. 13-23% of the FN and 22-26% of the present taxa are in-silico strains, so the model learns
  part of what a strain's spectrum looks like from a model that is off by 0.06-0.10. ONT shows no such gap: its
  errors swamp the spectrum (0.36-0.42 for strains at any identity).
- **Rescaling divergence by the gene factors does nothing**: per taxon, the raw mean divergence of the counted
  fragments separates FN from FP taxa at AUC 0.46-0.61, the mean of divergence over the gene's within- or
  between-species factor at 0.44-0.61.

### 9. The beyond-genus copies against the database's near pairs

[`incongruence.py`](incongruence.py) ([incongruence.txt](incongruence.txt)) checks section 7 from the database's
side: `gene_incongruence.tsv` (the v14 build's; v15 has the same reference genes) lists every two copies of a gene in
different genera within 0.05 k-mer distance, 105,390 copies with such a partner, 5,589 of them suspect. The copies the
counted fragments land on:

| fragments of | pe | se | pb | ont |
|---|---|---|---|---|
| FP, beyond the genus: copy has a copy of another genus within 0.05 | 0.17 | 0.17 | 0.22 | 0.21 |
| ... within 0.02 / 0.01 | 0.14 / 0.13 | 0.14 / 0.12 | 0.19 / 0.15 | 0.19 / 0.16 |
| ... and the read is no nearer the copy than that other genus's copy is | 0.12 | 0.12 | 0.14 | 0.19 |
| FP, same genus: copy has such a partner | 0.004 | 0.003 | 0.003 | 0.005 |
| own reads (FN and TP taxa): copy has such a partner | 0.002-0.003 | 0.002 | 0.002-0.007 | 0.003 |
| the FP taxon has 3 or more congeners in the training database | 0.75 | 0.77 | 0.70 | 0.70 |

- **A per-copy uniqueness scale would catch a fifth of this class.** 78-83% of the beyond-genus FP reads sit on
  copies that are unique in the database at 0.05, as section 7 found from the ZA tags. A rule "the read is no nearer
  the copy than another genus's copy is" removes 12-19% of these reads at a cost of 0.1-0.2% of own reads; with the
  class worth 10-15% of the FP at the knob (the "near novel" and "absent" rows of section 6's budget), that is about
  +0.0005 F1.
- **What could see the rest** is the copy's congruence with its own congeners: 70-77% of these FP taxa have three or
  more congeners in the database, so a copy whose distance to its nearest congener's copy is an outlier against the
  species' other genes (a transferred or contaminating copy sits far from the congeners' copies on that gene alone)
  can be flagged without any foreign copy. `ScanGene` computes that distance already and keeps only the verdict;
  writing both distances per copy (congener, foreign, with the rank) gives the graded scale and the congruence test in
  one pass.

### 10. A cross-check of section 6 with a plain converter ([genus_level.py](genus_level.py), [genus_level.txt](genus_level.txt))

A boosted classifier on the called taxa (all 108 features and the model's score; target "absent novel congener";
out-of-fold scores with samples grouped on the training rows, AUC 0.96-0.98, AP 0.60-0.72; the threshold of highest
species F1 there, applied once to the test rows). A converted absent taxon is no longer a false positive; a converted
present taxon is a false negative at species level.

| test rows | pe | se | pb | ont |
|---|---|---|---|---|
| TP / FP / FN, F1 at the knob | 26,932 / 835 / 1,610, 0.9566 | 4,271 / 134 / 202, 0.9622 | 2,048 / 39 / 51, 0.9785 | 21,514 / 1,023 / 1,176, 0.9514 |
| FP that are novel congeners | 718 (86%) | 99 (74%) | 33 (85%) | 894 (87%) |
| ceiling: all of them at genus level, F1 | 0.9689 | 0.9730 | 0.9863 | 0.9706 |
| converted: absent, genus right / wrong / present | 342 / 5 / 147 | 49 / 4 / 21 | 14 / 0 / 5 | 392 / 5 / 207 |
| F1 after conversion | 0.9598 (+0.003) | 0.9655 (+0.003) | 0.9806 (+0.002) | 0.9550 (+0.004) |
| "genus G is present" right / "unknown species of G" right | 0.99 / 0.69 | 0.93 / 0.66 | 1.0 / 0.74 | 0.99 / 0.65 |

A quarter to a third of the ceiling; the [2026-10-03 genus fallback](../2026-10-03-genus-fallback/README.md) found
+0.0002 on v5, and v15's model and features are the difference (the converter's top features:
`expected_gene_presence_unique_weighted`, `third_position_share`, `excess_high_share`, `gene_presence_ratio`,
`cluster_mean_ani`).

## What this means for classifying them

1. **New read features will not separate these errors.** The analogues computable here (congener ties and negative
   `ZA` values, failed seeds on related taxa, divergence spread over genes, gene-end and mate patterns) show at most
   AUC 0.6 between FP reads and the missed taxa's own reads, and nothing beyond divergence at the taxon level. That
   is a caution for the `fp-features` branch's similar features: they need an r226 build to show a gain.
2. **Strain alleles of wide species in the index.** FN taxa are strains, 89-97% of their fragments from a
   non-representative genome, of species with many genomes. Their reads would align at a closer allele, while a
   novel congener's would not. This attacks the identity overlap directly; the 2026-10-03 FN anatomy put its
   ceiling at about +0.009 pe F1.
3. **A per-species marker radius in place of `cluster_ani_radius`.** At build time, take the marker identity of each
   cluster's most distant genome to its representative, so the model can ask whether reads at 0.96 are within what
   the species' own strains show. `cluster_min_ani` is a genome-wide proxy (AUC 0.55); the radius as it is is
   constant.
4. **Genus-level calls for novel congeners (section 6), as an added output rather than a replacement.** The build
   would train a second model on the same rows (target: absent, in a genus with a novel congener) and ship it in the
   database. A run would then add "unknown species of G" lines for genera with reads and no confident call, and
   relabel the few calls the second model is sure of. On the simulated test sets: species F1 +0.002 to +0.003, genus
   calls 82-92% right; 64% where novel congeners are rare. Check how many such lines a real gut sample gets before
   making it a default; that is a reporting choice for you.
5. **The beyond-genus singletons (sections 7 and 9)** show nothing in the read, and little in the database: for
   98-99% of these reads no other genus's copy is within 2 edits. A graded uniqueness score over the database's copies
   would catch at most the 12-19% that section 9's rule removes. Three things could reach more:
   - a congruence test of each copy against its own congeners' copies (section 9);
   - a scan of whole genomes;
   - a per-copy foreign-read rate learned from the build's simulations (a fifth to two fifths of them).

   The class is 8-20% of FP taxa, so each is worth doing only if cheap.
6. **Test `--align_top` above 3 in crowded genera (section 5).** As many of a missed taxon's own reads land on a
   congener, at the same identity, as on the taxon itself, and the own species was never aligned against. The share
   grows with genus size. A cluster re-run of some soil test samples with `-c 10` would show whether a wider candidate
   list keeps these reads (and costs how much time), or whether the strains truly sit nearer the congener.
7. **ONT's alignment thresholds.** Only Nanopore loses own reads at the alignment itself: 15%, and those reads sit at
   0.90 identity where they do align. Which threshold rejects them (score, x-drop, identity) is not in the SAMs; it is
   worth a look before the next ONT model.

## Further approaches, ranked by what these data support

Beyond the items above (strain alleles, the marker radius, genus-level calls, `--align_top`, ONT's thresholds), with
the expected gains as test-set species F1 for pe unless said.

1. **Report novel congeners at genus level** (sections 6 and 10): +0.002 to +0.004 now, ceiling +0.008 to +0.019; a
   second model in the database and a `call_level` column (species / genus) in the output. The genus is right
   93-100% of the time; a third of the conversions demote a present species to its genus. Whether that trade reads as
   a gain to users is your call.
2. **Count the ambiguous strain reads** (section 5): 14-30% of FN taxa lose half or more of their aligned own reads
   to a close congener that is not called either; the pair splits the reads and the MAPQ filter drops them from
   both. The EM resolves most of them (`em_own_share` 1.0) and the `unfiltered` features (`fragments_all`,
   `em_fragments`) are in the model, but the count that drives the call is the filtered one. Two forms: EM-weighted
   fragments summed over the taxon and its congeners within 0.02 (`species_neighbours.tsv` has them) as a "complex
   depth" feature, so that the complex counts even where the species is ambiguous; and a joint call of such pairs
   ("*S. mitis* / *S. pneumoniae* group") where the reads do not separate them, instead of one miss and one absence.
   Testable first on the cluster with `--align_top` (item 6 above), which tells how much is candidate selection.
3. **Give the in-silico strains a real strain's spectrum** (section 8): fit `--omega` (and the per-frame rate gamma)
   so that the in-silico strains' `third_position_share` matches the real strains' per identity band (0.70 at
   0.95-0.97, 0.63 at 0.97-0.98 for pe reads); or draw omega per strain from the real strains' spread. One parameter
   in `insilico_strains.py`, no protal change; it removes a spectrum the model learns from 22-26% of its present
   taxa that no real organism shows. The F1 effect is unknown and probably small; the feature's meaning becomes right.
4. **Strain alleles for re-scoring, not for seeding** (a cheaper form of item 2 above). Store the other genomes'
   marker copies as differences against the representative (strains differ at ≤5% of positions; a few genomes per
   species chosen for diversity; about a gigabyte at r226) and re-align only the reads of taxa in the ambiguous band
   (identity 0.95-0.985) against them: a strain's reads gain (0.965 → 0.99 on the right allele), a novel congener's
   do not. Features: the identity gain on the best allele, the share of reads improved. The index and the seeding are
   untouched, the runtime cost is confined to band taxa. The ceiling of strains detected like representatives was
   +0.009; the FP side gains too, as the band empties of strains. The simulation must draw its strain genomes from
   genomes that are *not* in the allele table, or the gain is leakage.
5. **Placement on the terminal branch.** A strain attaches at the tip of its species' branch, a novel congener lower
   down; the sites where the species differs from its nearest congeners (its derived states) tell where: a strain's
   reads agree with the species there, a congener that branched at fraction *x* of the branch agrees at about *x* of
   them. protal's unique k-mers are a proxy of this already (`lu_per_kb` separates the band at AUC 0.80-0.88, which is
   why `ZU` showed nothing *more* at the boundary), but k-mer survival also falls with every sequencing error and with
   the organism's own mutations. Two cleaner forms, both needing the congeners' gene copies at build time (the build
   sketches them already): (a) per copy, the positions where it differs from its nearest congeners' copies, and at
   run time the share of those positions the taxon's reads agree with the species at; (b) "ancestral decoys": one
   extra reference per species placed along its branch toward its congeners (its derived states partly reverted), so
   that a novel congener's reads prefer the decoy and become "unknown species near T" while a strain's reads stay with
   T. Testable offline with these records plus the gene sequences of the FP/FN taxa and their three nearest
   neighbours from the training database (a cluster extraction; the reads' sequences are in the SAMs).
6. **A per-copy congruence table at build time** (section 9), in place of the binary suspect flag: per copy the
   distance to the nearest congener's copy and to the nearest copy of another genus (with its rank), from the pass
   that exists. A copy far from its congeners' on one gene alone is flagged without any foreign copy, which reaches
   the 80% of beyond-genus FP copies that are unique in the database; the graded distances become read weights or
   taxon features. The class is worth 10-15% of the FP.
7. **Priors from real data.** Most FN and half of the FP have 1-2 fragments, where no read feature can help and the
   sample's depth already does what a prior can in simulation. A per-species prevalence by habitat, from real
   profiles (a habitat inferred from the sample's confident calls, then the species' prior in that habitat as a
   feature), is the one source of information the simulation cannot supply and cannot validate; it needs real
   benchmarks (mock communities, CAMI, your own collections) and is where the largest residual sits.

Not worth pursuing, as measured: more read-level features of the alignments at the boundary (AUC 0.57-0.64); the
spread or gene-factor rescaling of divergence (sections 4 and 8); a database-internal uniqueness scale on its own for
the beyond-genus class (sections 7 and 9).

## Implemented and prepared from the list

- **In-silico strain spectrum (item 3): done.** `scripts/insilico_strains.py --omega auto` (the default now) samples
  200 of the table's real strains, compares their marker-gene copies with their representatives' along the shared
  12-mers (`substitutions`, no alignment; a stretch past an indel is skipped), measures the share of substitutions
  on third codon positions, and takes the omega at which `mutate()` reproduces that share on the representatives'
  genes (`calibrate_omega`, a grid of 0.01-1.0 interpolated in log omega). The build's `insilico_strains.log` says
  what was measured and chosen; `--omega 0.15` keeps the old behaviour. Tests: `scripts/test_insilico_strains.py`
  (9, with the substitution finder and the calibration). `docs/databases.md` updated. The next build's in-silico
  strains will carry a real strain's spectrum; whether F1 moves is for that build to show.
- **Placement on the terminal branch and the species' alleles (items 4 and 5, merged): measured by the build.**
  `scripts/ancestry_sites.py` (`scripts/test_ancestry_sites.py`, 3 passing; first written here as
  `derived_sites.py`, whose cluster run did not work, so the build now runs it where the inputs are:
  `build_gtdb_database.py` with `--share-logs`, after the error reads, into `model_logs/ancestry_sites/<read
  type>.{summary.txt,auc.tsv,taxa.tsv.gz,fragments.tsv.gz}`, the training database's full reference kept until
  then and its genes unpacked for it; the pipeline test `test_f_scenarios` covers it) takes the error SAMs of
  one read type, the training database's `reference.fna`, its build input `full_reference.fna.zst` (every genome's
  marker genes under the species' taxid), the taxonomy and the held-out list. For every FP and FN taxon it finds the
  sites where its gene copies differ from its nearest congener (and from the majority of the three nearest) and the
  sites where the species' own alleles differ from the representative (polymorphic), and classes the congener
  sites as *fixed* in the species (every known allele shares T's base: a synapomorphy) or polymorphic. For every
  counted read it counts the covered sites of each class and the read's base there. Signals per taxon: the share
  of the species' base at the congener sites and at the fixed sites (a strain near 1, a novel congener at the
  fraction of T's branch it shares), the share of mismatches at polymorphic sites and the *fixed-site identity*
  (identity over the non-polymorphic bases: a strain's mismatches fall where the species varies, a congener's do
  not). Output: per record, per taxon, and the AUC of each signal for the missed taxa's own reads against the FP
  reads, overall and within identity bands. By hand, where the gene sequences are (the full reference is one pass
  of 86 GB, about 10 minutes):

      protal --db training_db/database.protal --unpack_db --unpack_dir /tmp/training_db_files
      python3 scripts/ancestry_sites.py --sams OUTDIR/model_logs/error_reads/pe --reference /tmp/training_db_files/reference.fna \
          --full-reference training_db/full_reference.fna.zst --taxonomy OUTDIR/internal_taxonomy.dmp \
          --heldout OUTDIR/heldout_species.txt --out ancestry_pe

  The next build's `model_logs/ancestry_sites/pe.summary.txt` says it: an AUC of the fixed-site signals well above
  identity's within the 0.96-0.985 bands is the go for the alleles' half in protal (below); the congener sites are
  in protal already, as features the same build's models carry.
- **Strain alleles for re-scoring (item 4): merged into the above.** No download is needed: the build already has
  every GTDB genome's marker genes. Re-aligning reads to alleles collapses into the polymorphic-site mask: a mismatch
  at a site where the species' own strains vary is forgiven, one at a fixed site counts, and at the fixed sites that
  separate the species from its congeners the read's base says which side it is on. The database stores per copy
  only the sparse polymorphic-site mask (species with two or more genomes; a few hundred MB at r226, read per
  species on demand); the congener sites are computed at run time from the gene store and `species_neighbours.tsv`
  for the genes a taxon's reads hit.
- **In protal, first half (0.7.9, 2026-10-07): the congener sites, no database change.** `src/SequenceUtils/AncestrySites.h`
  pairs a species' gene copy with its nearest congener's along their shared 12-mers (the same method as
  `ancestry_sites.py`), keeps the differing positions with the congener's base (a few hundred bytes), computes them
  once per run for each (species, gene) a record touches (`GenomeLoader::AncestrySitesOf`, a thread-safe cache capped
  at 200,000 pairs) and counts per best record, from its CIGAR, the sites it covers, those where the read has the
  species' base (M runs) and those where it has the congener's (X runs, the read's base looked up). Three features
  in a new default group `ancestry` (`model_features.py`, [features.md](../../features.md)):
  `ancestry_sites_per_record`, `ancestry_agreement`, `ancestry_congener_share` (-1 without a site). Unit tests in
  `tests/test_AncestrySites.cpp` (4; the suite 441 of 442, the failure `IlluminaSimulation.AFailedStreamIsCutOff`
  being in HEAD before these changes), the 142 end-to-end tests (which check that every default feature is in the
  dump) and the trainer's 41 feature-set tests pass. On the mini database's mock community (20,000 pairs): *Mockella
  alpha* at identity 0.988 covers 10.8 sites per record and agrees with its reference at 0.957 of them (congener's
  base 0.036), *Mockella beta* at 0.996 agrees at 0.998, *Fakibacter gamma*, without a congener in the database,
  gets -1. Runtime: the comparison runs for the genes a taxon's reads hit only (tens of microseconds each); per
  record a binary search per CIGAR run, the sites taken from a per-chunk cache in front of the run's. The polymorphic
  mask (the alleles' half) waits for the cluster measurement to say whether the fixed-site refinement pays.

## Also here, not run

[`error_read_features.py`](error_read_features.py) (with [`run_cluster.sh`](run_cluster.sh) and
[`test_error_read_features.py`](test_error_read_features.py), 4 tests passing) would take the FP and FN reads again
from a build's samples on its scratch disk. It also writes read features for every row, TP and TN included.
[`analyse_error_reads.py`](analyse_error_reads.py) analyses that output. They were not needed here, because the logs
held the SAMs. They are what a fair FP-against-TP comparison at the read level would need, instead of the biased
"own (TP)" group above.
