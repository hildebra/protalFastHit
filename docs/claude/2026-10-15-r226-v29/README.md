# The r226 v29 build: the gates hold on deep samples, GBM beats the default forest, the gate models cost an hour

Date: 2026-10-15. Branch `audit-fixes` at `55f850d`. Nothing in protal or the scripts was changed for this report.

**Data.**
- **v29:** SLURM 24261409 on q1536n3 (1.5 TB), 2026-10-15 13:29 to 17:45, 120 threads, `--scratch` on the node's
  SSD. protal, the simulator and the scripts at `4d89d86` (0.8.4 with the build audit's fixes `3305b2c`;
  `protal_db/build_metadata.tsv`). The command: `build_gtdb_database.py --inputs .../refGenos/ ... --share-logs
  --test-c-gates -t 120 --scratch ...`. Everything else is a default of `4d89d86`:
  - database gates `pe:4,se:4,pb:16,ont:8`, picked by stride;
  - models also at `pe:1/8,se:1/8,pb:8/32,ont:1/16` (`--model-gates`);
  - forest models as the database's default, GBM models beside them (`--model forest,gbm`);
  - paired-end training at 75, 100, 150 and 250 bp; single-end Illumina at 150 and 250 bp, and Ultima reads;
  - the 98-feature default set.
- **The archive:** `local/protal0.8.4_r226_v29_share.tar.gz`, unpacked to `local/v29/` (not in git). It holds the
  console, `logs/`, `model_logs/` and the database gate's training and test tables (`work/`).
- **Compared with v27** ([2026-10-13-r226-v27-c-gates](../2026-10-13-r226-v27-c-gates.md); unpacked in `local/v27/`):
  - Built at `9c27756` (0.8.3) on the same node: gates `pe:1,se:1,pb:8,ont:1` by hash, GBM models only, 150/250 bp
    paired-end training.
  - Same download, seed, genome table, hold-outs and scenarios as v29. The scenario samples are the same
    simulations. The design samples are not, because the read setups changed.
  - v28, the run before v29, crashed at the parity check (fixed in `b42cbc3`; [build audit](../2026-10-15-gtdb-build-audit/README.md)).
- **Commands** (WSL, `taskset -c 0-1 nice -n 5`):
  - [pairs.py](pairs.py) `local/v29/.../model_logs local/v27/.../model_logs` → [pairs.txt](pairs.txt). It reads the
    trainers' `*.calls.tsv.gz` (every row's score, knob and call as protal calls). It gives:
    - forest against GBM on the same rows, with a bootstrap over samples and the true positives at matched false
      positives;
    - v29 against v27 on the test samples both builds simulated alike;
    - the paired-end test samples by read length.
  - awk one-liners over `c_gates*.tsv`, the reports' knob and marker-truth sections, `logs/cpu_jobs.tsv` and
    `logs/protal_cpu_{training,test}.tsv` of both builds (sections 3-6).

## Summary

- **The gates hold on deep, complex samples.** v27's report left this open: no streamed scenario sample had been
  profiled at a higher gate. On the 20 scenario hold-out samples per read type that both builds simulated alike,
  v29's GBM at its gates (pe C 4, pb C 16, ont C 8, by stride) is level with v27's GBM at C 1 / C 8 by hash:
  - pe +0.0007, 95% interval [−0.0003, +0.0019];
  - pb +0.0013 [+0.0003, +0.0024];
  - ont +0.0013 [+0.0006, +0.0021].

  Among the test samples' gate rows (without `--flex_indels`), the database's gate is the best or within 0.0003 of the
  best row for pe, pb and ont. It saves 23% (pe) and 39% (ont) of the alignment CPU against C 1, and 49% (pb) against
  no gate. The models' other gates also profiled the streamed scenario samples. On those samples, every gate from C 1
  to C 16 (pb C 32), with its own GBM model, is within 0.003 of the database's gate. The exception is pe C 8, whose
  model calls at knob 0.5.
- **The database's default model kind, the forest, is the weaker one.** The GBM models beat the forest in every read
  type and set, at each kind's own knob:
  - species held out +0.003..+0.009;
  - scenario hold-out +0.005..+0.009;
  - test design +0.001..+0.004.

  At equal false positives, GBM finds more true positives in 58 of 60 comparisons. The exceptions are the paired-end
  gut samples and five pb test taxa. GBM's extra gut misses are 41 of 43 in one sample, mostly species with 1-2 marker
  reads. Against v27's GBM, v29's forest is 0.0045-0.0074 lower on the shared scenario samples.

  0.8.4 made the forest the default (`ce2940e`) on small-world and real-mock evidence. The real mocks found the
  kinds equal in F1 ([mock sparse gates](../2026-10-14-mock-sparse-gates/README.md)). At r226, GBM is better.
- **The models at other gates cost about 70 minutes and gain nothing measurable as models.** Their held-out F1 is
  within 0.0009 of the database gate's model for each read type and kind. Their rows are useful as a test: they are
  the only rows of the streamed deep samples at other gates. Profiling both collections at two more gates per
  read type made up 48% of the collection's protal CPU, against 21% for the database gate's own rows. The
  collections took 2:27 against v27's 1:06; the whole run took 4:16 against 2:38.
- **The trainer's knob rule reset 7 of the 24 models to knob 0.5:**
  - pe C 8 GBM;
  - se GBM at all three gates;
  - se C 8 forest;
  - ont GBM at C 8 and C 1.

  In each, the bootstrap's knob won in 96-100% of the resamples. The rule dropped it for losing 0.0001-0.0018 of F1 on
  the weighted test set (audit finding E2). The test set's design samples always want a lower knob than the held-out
  species. For pe C 8 GBM, 0.78 would gain 0.0034 on the scenario hold-out and lose 0.0038 on the design samples. At
  0.5 that model calls 183 false positives per soil_shallow sample. `f1-integration` no longer lets the test set drop
  the knob.
- **Unchanged since v27:**
  - Present species with marker reads but no row: ONT 3.7% (half of ONT's marker-truth misses), pe 1.6%, pb 1.4%.
  - PacBio's composition: species' depth +27.6% (v27 +21.0%), unknown share −0.18.
  - Parity: every model, the gate models included, gives protal's probabilities and features.
- **The grid's sparse index (emulated):**
  - F1: within 0.0005 of the stride at pe C 4 and pb C 8, but −0.0020 at pb C 16 and −0.0046 at pe C 8.
  - CPU: the console's figures are those of the emulation on the whole index (pe C 4: 2.8× the stride's), not of a
    sparse index.

## 1. Headline numbers against v27

F1 from the consoles. Species held out at knob 0.5; the test set and the soil hold-out samples as protal calls (each
model at its own knob). v27's models are GBM, its default then.

| read type | | species held out | test set | soil hold-out | soil_shallow hold-out |
|---|---|---:|---:|---:|---:|
| pe | v27 GBM (C 1) | 0.9490 | 0.9651 | 0.9529 | 0.9312 |
| | v29 GBM (C 4) | 0.9491 | 0.9647 | 0.9540 | 0.9311 |
| | v29 forest (C 4, the default) | 0.9406 | 0.9634 | 0.9450 | 0.9250 |
| se | v27 GBM (C 1) | 0.9572 | 0.9634 | 0.9525 | — |
| | v29 GBM (C 4) | 0.9573 | 0.9668 | 0.9531 | — |
| | v29 forest (C 4) | 0.9487 | 0.9625 | 0.9455 | — |
| pb | v27 GBM (C 8, hash) | 0.9712 | 0.9835 | 0.9684 | 0.9609 |
| | v29 GBM (C 16) | 0.9718 | 0.9847 | 0.9703 | 0.9618 |
| | v29 forest (C 16) | 0.9670 | 0.9819 | 0.9604 | 0.9573 |
| ont | v27 GBM (C 1) | 0.9643 | 0.9729 | 0.9588 | 0.9493 |
| | v29 GBM (C 8) | 0.9656 | 0.9751 | 0.9604 | 0.9517 |
| | v29 forest (C 8) | 0.9588 | 0.9735 | 0.9469 | 0.9448 |

The species-held-out and test columns compare different samples: v29 trains and tests paired-end reads of 75-250 bp,
v27 150 and 250 bp; se's design changed too (261 training samples against 339). Only the scenario samples are the
same simulations in both builds (section 3).

## 2. Forest against GBM

Paired on the same rows of v29 ([pairs.txt](pairs.txt), section 1). Each kind is scored at its own knob, as protal
calls (pe forest 0.69, GBM 0.74; se 0.58 and 0.5; pb and ont 0.5 both). The interval is 95% over 1,000 resamples of
the samples.

| read type | set | samples | F1 forest | F1 GBM | GBM − forest | 95% interval | TP at GBM's FP: forest/GBM | TP at forest's FP: forest/GBM |
|---|---|---:|---:|---:|---:|---|---|---|
| pe | held out, design | 308 | 0.9644 | 0.9697 | +0.0053 | [+0.0039, +0.0068] | 15,036 / 15,259 | 15,217 / 15,374 |
| pe | held out, scenarios | 50 | 0.9420 | 0.9513 | +0.0094 | [+0.0080, +0.0110] | 56,610 / 57,736 | 57,236 / 58,130 |
| pe | test, design | 104 | 0.9634 | 0.9647 | +0.0013 | [−0.0016, +0.0042] | 5,536 / 5,590 | 5,625 / 5,666 |
| pe | test, scenarios | 20 | 0.9477 | 0.9530 | +0.0053 | [+0.0015, +0.0079] | 22,464 / 22,761 | 22,849 / 23,051 |
| pe | test, gut | 4 | 0.9851 | 0.9766 | −0.0085 | [−0.0312, +0.0057] | 1,475 / 1,479 | 1,516 / 1,498 |
| se | held out, design | 231 | 0.9632 | 0.9679 | +0.0047 | [+0.0031, +0.0064] | 11,405 / 11,519 | 11,416 / 11,522 |
| se | held out, scenarios | 30 | 0.9459 | 0.9541 | +0.0082 | [+0.0064, +0.0101] | 36,014 / 36,621 | 35,920 / 36,542 |
| se | test, design | 78 | 0.9625 | 0.9668 | +0.0042 | [+0.0011, +0.0072] | 4,059 / 4,100 | 4,097 / 4,141 |
| se | test, scenarios | 12 | 0.9542 | 0.9599 | +0.0058 | [+0.0029, +0.0079] | 14,189 / 14,357 | 14,281 / 14,437 |
| pb | held out, design | 186 | 0.9819 | 0.9849 | +0.0030 | [+0.0014, +0.0046] | 6,024 / 6,049 | 6,011 / 6,047 |
| pb | held out, scenarios | 50 | 0.9650 | 0.9700 | +0.0050 | [+0.0039, +0.0059] | 43,201 / 43,714 | 43,474 / 43,909 |
| pb | test, design | 40 | 0.9819 | 0.9847 | +0.0027 | [+0.0005, +0.0051] | 1,471 / 1,477 | 1,467 / 1,462 |
| pb | test, scenarios | 20 | 0.9677 | 0.9738 | +0.0061 | [+0.0041, +0.0080] | 17,944 / 18,241 | 18,173 / 18,359 |
| ont | held out, design | 186 | 0.9753 | 0.9803 | +0.0050 | [+0.0031, +0.0068] | 6,102 / 6,187 | 6,144 / 6,206 |
| ont | held out, scenarios | 50 | 0.9566 | 0.9636 | +0.0070 | [+0.0060, +0.0081] | 43,490 / 44,151 | 43,746 / 44,319 |
| ont | test, design | 41 | 0.9735 | 0.9751 | +0.0016 | [−0.0015, +0.0047] | 1,485 / 1,489 | 1,486 / 1,489 |
| ont | test, scenarios | 20 | 0.9574 | 0.9661 | +0.0087 | [+0.0057, +0.0116] | 18,068 / 18,399 | 18,252 / 18,521 |

- **GBM ranks the taxa better, not only at a better knob.** pairs.txt has 30 sets: the table's and each scenario
  alone. In each set, each kind gets the other's false-positive count. In 58 of these 60 comparisons GBM then finds
  more true positives. The two exceptions:
  - pb's test design rows: 1,462 against 1,467 at the forest's 15 FP;
  - pe gut (below).

  The largest gaps are on the soil samples: pe +0.0090, pb +0.0099, ont +0.0135. There the forest calls 1.7-2.2×
  GBM's false positives for no more true positives (pe +62, pb −20, ont −40).
- **The exception, pe gut (4 hold-out samples):** GBM misses 66 present species, the forest 29.
  - 41 of GBM's 43 extra misses are in one sample (`sc_gut_pe_p20000000_s_4`), 34 of them with 1-2 marker reads. The
    forest calls them at 0.8-0.9, GBM at 0.1-0.5.
  - The one abundant species among them is *Undibacterium parvum* (689 marker reads, ~1× depth) in sample 1: GBM
    scores it 0.65 against its knob of 0.74.
  - The real D6331 mock showed the same pattern (GBM missed two abundant paired-end species). Too few samples to call
    it a gut weakness of GBM, but worth watching on the real mocks.
- **The forest's held-out deficit is not new.** v27's GBM matched v29's GBM. The forest is the default because
  `ce2940e` chose it on the small worlds and the strain world. At r226 every read type prefers GBM.

## 3. The gates on the same scenario samples, and on the test samples' gate rows

### v29 against v27 on the scenario hold-out samples

The test samples whose present taxa overlap by 95% or more between the builds ([pairs.txt](pairs.txt), section 2). A
present taxon without a row in one build counts as a false negative there.

| read type | samples | present | rowless v27/v29 | F1 v27 GBM | F1 v29 GBM | v29 GBM − v27 | 95% interval | F1 v29 forest | v29 forest − v27 |
|---|---:|---:|---|---:|---:|---:|---|---:|---:|
| pe | 20 | 24,247 | 0 / 1 | 0.9522 | 0.9529 | +0.0007 | [−0.0003, +0.0019] | 0.9477 | −0.0045 |
| se | 8 | 10,733 | 59 / 60 | 0.9629 | 0.9654 | +0.0025 | [+0.0006, +0.0058] | 0.9584 | −0.0045 |
| pb | 20 | 18,882 | 4 / 0 | 0.9725 | 0.9738 | +0.0013 | [+0.0003, +0.0024] | 0.9677 | −0.0048 |
| ont | 20 | 19,179 | 0 / 9 | 0.9645 | 0.9658 | +0.0013 | [+0.0006, +0.0021] | 0.9572 | −0.0074 |

- The scenario samples are identical simulations (same present taxa). They cover gut, moderate, soil, soil_shallow and
  host at 5-20 M pairs or 1.5-6 Gb. On them, v29's GBM is level with v27's GBM or slightly better. It runs pe at C 4
  instead of C 1 and ONT at C 8 instead of C 1, has new features and was trained on more read lengths.
- The gates make few species rowless. pe loses 1 of 24,247 at C 4 and ONT 9 of 19,179 at C 8. pb at C 16 by stride
  loses none, where v27's C 8 by hash lost 4.
- The design samples cannot be compared: the read setups changed, and 0 of the 77 shared pe sample names have the
  same taxa.

### The scenario hold-out samples at the models' gates

The models' other gates profiled every sample again, the streamed deep ones too, and each gate's model scored its
own rows. This is the only look at the deep samples off the database's gate. GBM models, F1 as protal calls (from
`summary.txt`):

| read type | gate (knob) | gut | moderate | soil | soil_shallow |
|---|---|---:|---:|---:|---:|
| pe | 1 (0.75) | 0.9739 | 0.9706 | 0.9535 | 0.9284 |
| pe | **4 (0.74)** | 0.9766 | 0.9699 | 0.9540 | 0.9311 |
| pe | 8 (0.5) | 0.9791 | 0.9645 | 0.9502 | 0.9234 |
| se | 1 (0.5) | — | 0.9657 | 0.9544 | — |
| se | **4 (0.5)** | — | 0.9673 | 0.9531 | — |
| se | 8 (0.5) | — | 0.9685 | 0.9558 | — |
| pb | 8 (0.5) | 0.9917 | 0.9821 | 0.9704 | 0.9617 |
| pb | **16 (0.5)** | 0.9924 | 0.9828 | 0.9703 | 0.9618 |
| pb | 32 (0.5) | 0.9913 | 0.9819 | 0.9695 | 0.9610 |
| ont | 1 (0.5) | 0.9888 | 0.9774 | 0.9597 | 0.9503 |
| ont | **8 (0.5)** | 0.9906 | 0.9782 | 0.9604 | 0.9517 |
| ont | 16 (0.5) | 0.9895 | 0.9775 | 0.9611 | 0.9528 |

- The gate changes little on the deep samples: all within 0.003 of the database's gate, except pe C 8, where the
  model's knob of 0.5 adds false positives (section 5). One gut species is ~0.0007 of F1 for pe.
- se C 8 and ont C 16 hold on them: se moderate +0.0012 and soil +0.0027; ont soil +0.0007 and soil_shallow +0.0011.

### The test samples' gate rows (`c_gates_gbm.tsv`)

The non-streamed test samples are profiled again at each gate in the same protal run. Each row is scored by the model
of its nearest gate (`model`). Since `3305b2c` a present species with marker reads but no row counts as a false
negative (`unseen`). Δ is the F1 difference from the database's gate. CPU is the alignment's CPU seconds summed over
the samples.

| read type | gate | pick | model | unseen | TP | FP | FN | F1 | Δ | CPU s |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| pe | off | | 1 | 367 | 19,266 | 631 | 1,648 | 0.9442 | −0.0014 | 3,575 |
| pe | 1 | | 1 | 367 | 19,303 | 645 | 1,611 | 0.9448 | −0.0008 | 2,681 |
| pe | 2 | stride | 1 | 367 | 19,306 | 650 | 1,608 | 0.9448 | −0.0008 | 2,046 |
| pe | **4** | **stride** | 4 | 369 | 19,383 | 701 | 1,531 | **0.9456** | | **2,053** |
| pe | 4 | hash | 4 | 399 | 19,384 | 726 | 1,523 | 0.9452 | −0.0004 | 1,840 |
| pe | 4 | grid | 4 | 362 | 19,377 | 706 | 1,544 | 0.9451 | −0.0004 | (5,799) |
| pe | 8 | stride | 8 | 372 | 19,700 | 1,246 | 1,207 | 0.9414 | −0.0042 | 1,748 |
| pe | 8 | grid | 8 | 420 | 19,706 | 1,444 | 1,216 | 0.9368 | −0.0088 | (5,092) |
| pe | 16 | stride | 8 | 403 | 19,703 | 1,228 | 1,196 | 0.9421 | −0.0035 | 1,669 |
| se | 1 | | 1 | 235 | 18,510 | 687 | 1,026 | 0.9558 | +0.0003 | 3,183 |
| se | **4** | **stride** | 4 | 241 | 18,457 | 644 | 1,078 | **0.9554** | | **2,325** |
| se | 8 | stride | 8 | 254 | 18,547 | 693 | 988 | 0.9566 | +0.0012 | 1,931 |
| se | 16 | stride | 8 | 296 | 18,560 | 829 | 968 | 0.9538 | −0.0016 | 1,888 |
| pb | off | | 8 | 121 | 5,899 | 129 | 395 | 0.9575 | −0.0003 | 2,404 |
| pb | 8 | stride | 8 | 122 | 5,903 | 131 | 391 | 0.9577 | −0.0001 | 1,471 |
| pb | **16** | **stride** | 16 | 123 | 5,909 | 136 | 385 | **0.9578** | | **1,232** |
| pb | 16 | hash | 16 | 133 | 5,899 | 135 | 395 | 0.9570 | −0.0008 | 1,238 |
| pb | 16 | grid | 16 | 131 | 5,915 | 158 | 390 | 0.9557 | −0.0020 | (1,870) |
| pb | 32 | stride | 32 | 128 | 5,884 | 125 | 410 | 0.9565 | −0.0013 | 1,075 |
| ont | 1 | | 1 | 361 | 6,063 | 228 | 682 | 0.9302 | −0.0003 | 2,005 |
| ont | **8** | **stride** | 8 | 365 | 6,062 | 222 | 683 | **0.9305** | | **1,218** |
| ont | 8 | hash | 8 | 378 | 6,051 | 218 | 693 | 0.9300 | −0.0005 | 1,222 |
| ont | 16 | stride | 16 | 371 | 6,061 | 219 | 682 | 0.9308 | +0.0003 | 1,087 |
| ont | 16 | grid | 16 | 415 | 6,063 | 296 | 665 | 0.9266 | −0.0040 | (1,748) |
| ont | 32 | stride | 16 | 420 | 6,014 | 226 | 726 | 0.9267 | −0.0039 | 935 |

The forest's table (`c_gates.tsv`) agrees on the database gates and the picks. These F1s are lower than the test set's
because the rows pool the design samples with the non-streamed scenario samples (soil_shallow and others) and count
the unseen species.

- **The database gates are right.**
  - pe C 4 is the best paired-end row. It aligns with 23% less CPU than C 1 and 43% less than no gate.
  - pb C 16 and ont C 8 are within 0.0003 of the best row. Against no gate or C 1 they save 49% (pb) and 39% (ont)
    of the CPU.
  - Beyond them, ont C 16 is free (−11% CPU) and se C 8 is free (−17%). pb C 32 (−0.0013) and ont C 32 (−0.0039)
    start to lose.
- **pe C 8 and C 16 rows reflect the C 8 model's knob, not the gate** (section 5): +545 FP against C 4, and 324 FN
  fewer. The forest's C 8 model (knob 0.67) scores pe C 8 at 0.9426 against 0.9419 at C 4.
- **Stride against hash, at the database gates:** the hash loses 7-30 more present species (unseen 399 against 369
  for pe, se 248 against 241, pb 133 against 123, ont 378 against 365) and 0.0004-0.0008 of F1. This confirms v27's
  choice.
- **The grid's sparse index (emulated):**
  - It costs ≤ 0.0005 at pe C 4 and pb C 8, −0.0020 at pb C 16, −0.0042 at ont C 16 and −0.0046 at pe C 8 against
    the stride with the same model.
  - It adds rows (pe C 8: 96,905 taxa against 89,954) and false positives (+198 at pe C 8). That fits the extra
    candidates of single-seed anchors.
  - These models were trained on stride rows, so a model trained on grid rows might recover some of this.
  - The CPU figures in parentheses are the emulation's: the whole index is loaded and every core looked up
    (`--sparse_index emulate`; `build_gtdb_database.py --help` says so). They say nothing about a sparse index's
    speed, but the console prints them beside the others without that note.
- **`--flex_indels`:**
  - ONT: +0.0007..+0.0008 at C 4-8, for 14-20% more CPU.
  - se's Ultima samples: 0.9559-0.9567 at C 0-4, against 0.9557 without it at the database's gate (the console's
    line). That is +0.0002..+0.0010; the forest's rows show −0.0001..−0.0002. Its CPU cannot be separated from the
    gate's: these samples have no row without it except at the database's gate.

  Neither is worth turning on by default.

## 4. The models at other gates

| read type | kind | species held out F1 at the gates (database's in bold) |
|---|---|---|
| pe | forest | C 1 0.9399, **C 4 0.9406**, C 8 0.9408 |
| pe | GBM | C 1 0.9492, **C 4 0.9491**, C 8 0.9496 |
| se | forest | C 1 0.9486, **C 4 0.9487**, C 8 0.9487 |
| se | GBM | C 1 0.9564, **C 4 0.9573**, C 8 0.9569 |
| pb | forest | C 8 0.9672, **C 16 0.9670**, C 32 0.9673 |
| pb | GBM | C 8 0.9720, **C 16 0.9718**, C 32 0.9724 |
| ont | forest | C 1 0.9586, **C 8 0.9588**, C 16 0.9589 |
| ont | GBM | C 1 0.9657, **C 8 0.9656**, C 16 0.9663 |

- Each model is trained and held out on its own gate's rows. All are within 0.0009 of the database gate's model. The
  features do not shift with the gate in any way the models can use; the `length_free` features and the unique
  k-mers counted from the alignment were built for that. v27 already found the C 1 model scoring C 2-4 rows without
  loss.
- This build cannot score a C 1 row with the C 4 model: the gate tables are not in the archive. So the direct test of
  what a gate model adds is not possible here. Nothing in the numbers suggests it adds anything.
- **What they cost.** The alignment's CPU seconds in both collections (`protal_cpu_*.tsv`):

  | rows | training | test | share |
  |---|---:|---:|---:|
  | the database's gate | 41,013 | 16,359 | 21% |
  | the models' other gates (`--model-gates`) | 94,950 | 37,127 | 48% |
  | the test-only gate rows (`--test-c-gates`: off, C 2, other picks, grid, flex) | — | 84,000 | 31% |
  | all | 135,963 | 137,486 | 273,449 |

  v27 aligned 110,000 CPU s in all (database gate 75,000, test gate rows 35,000). The database gate's rows got
  cheaper: −24%, on about as many samples (1,429 against 1,482). The two extra profilings per collection cost
  132,000 CPU s, and the collections took 2:27 against 1:06 (section 6). The grid's emulated rows are 34,000 CPU s of
  the test-only rows.

## 5. The knob rule

`choose_global_knob` takes the median of the bootstrap resamples' best thresholds if it beats 0.5 in 95% of them.
It then drops that knob if it loses any F1 on the weighted test set (the design's test samples plus the scenario
hold-out samples, the latter weighted 0.25). In v29 that rule reset 7 of the 24 models to 0.5. The other models
either kept their knob (pe forest and GBM at C 1 and C 4, pe C 8 forest, se forest at C 1 and C 4, pb C 32 forest) or
had no knob with 95% support (the rest of pb and ont). The Δ columns are F1 at the dropped knob minus F1 at 0.5:

| model | bootstrap knob (support) | test set F1 at it / at 0.5 | held out: Δ | test design: Δ | test scenarios: Δ |
|---|---|---|---:|---:|---:|
| pe C 8 GBM | 0.78 (100%) | 0.9569 / 0.9570 | +0.0058 | −0.0038 | +0.0034 |
| se GBM | 0.64 (96%) | 0.9630 / 0.9636 | +0.0010 | −0.0009 | −0.0002 |
| se C 1 GBM | 0.67 (100%) | 0.9623 / 0.9641 | +0.0019 | −0.0039 | +0.0004 |
| se C 8 GBM | 0.67 (97%) | 0.9639 / 0.9645 | +0.0016 | −0.0018 | +0.0007 |
| se C 8 forest | 0.60 (98%) | 0.9586 / 0.9587 | +0.0014 | −0.0017 | +0.0017 |
| ont GBM | 0.57 (97%) | 0.9677 / 0.9683 | +0.0007 | −0.0007 | −0.0005 |
| ont C 1 GBM | 0.52 (98%) | 0.9677 / 0.9679 | +0.0002 | −0.0013 | +0.0001 |

- Every dropped knob was higher than 0.5. Every one gained on the held-out species and lost on the design's test
  samples (−0.0007..−0.0039). On the scenario hold-out samples the effect is mixed (−0.0005..+0.0034). The test set
  wins because the design samples have most of its weight.
- The margins on the weighted test set are 0.0001-0.0018. For most of these models it hardly matters which knob is
  kept. pe C 8 GBM is the exception: it calls at 0.5 where its neighbours call at 0.74-0.75, with 183 false positives
  per soil_shallow sample against 102.
- The direction is systematic, not noise. The design's test samples prefer a lower knob than the species held out.
  Their species are in the final model's training samples, and the held-out species are not. v27 saw the same: its
  test set's best threshold was 0.51 against the knob's 0.73.
- The F1 audit's E2 is this rule. On `f1-integration` the test set is scored at the knob and at 0.5 for the report
  only; it no longer drops the knob. With it, these seven models would have kept the knobs above.
  `f1-integration` also changes the held-out scheme (`--cv-scheme`), so its knobs may differ.

## 6. Build time

Same node as v27 (q1536n3).

| phase | v27 | v29 |
|---|---|---|
| converting the release (background) | 0:29:58 | 0:46:11 |
| start to the collections (step 6/11) | 1:19:36 | 1:24:00 |
| both collections | 1:06:29 | 2:26:49 |
| `protal_db` (background, idle class) | 1:05:12 | 1:19:17 |
| training | 0:03:03 (GBM) | 0:08:38 (forest 4:30, then GBM 4:08; 12 models each) |
| run, to the last report | 2:38:20 | 4:12:36 |
| share archive | — | 0:03:36 (1.0 GB) |

- The collections are the difference: +80 minutes. The models' other gates are half of the collections' protal CPU,
  so roughly 70 of those minutes (section 4).
- The conversion read the release from `/hpc-home` 16 minutes longer than in v27 (network file system). It runs beside
  the genome table, so it delayed the collections by only ~4 minutes.
- Training both kinds at three gates takes 8.6 minutes. The parity of all 24 models takes 40 seconds.

## 7. Other observations

- **Present species without a row** (marker truth's `unseen`: marker reads, no record), share of the present species:

  | read type | v27 held out | v29 held out | v29 test set |
  |---|---:|---:|---:|
  | pe | 1.51% | 1.57% | 2.79% |
  | se | 1.61% | 1.57% | 2.72% |
  | pb | 1.38% | 1.37% | 1.26% |
  | ont | 3.63% | 3.72% | 5.16% |

  The stride gates did not raise these numbers: they were the same in v27 (C 1, pb C 8 by hash). For ONT they are
  1,946 of its 3,887 marker-truth misses with species held out. ONT recall has its largest reserve here, and the cause
  has not been looked at (the error-read SAMs hold only the rows' reads).
- **Read length** (pe test set, 26 design samples each, as protal calls; [pairs.txt](pairs.txt) section 3):

  | read length | F1 forest | F1 GBM | FP/sample forest | FP/sample GBM |
  |---:|---:|---:|---:|---:|
  | 75 | 0.9558 | 0.9615 | 1.81 | 1.19 |
  | 100 | 0.9605 | 0.9585 | 1.62 | 1.27 |
  | 150 | 0.9752 | 0.9747 | 1.00 | 0.81 |
  | 250 | 0.9642 | 0.9660 | 2.54 | 1.65 |

  Trained on 75 bp, the 75 bp samples call 1.5× the false positives of 150 bp, not the 6× of a 150 bp model
  ([length-free features](../2026-10-14-length-free-features/README.md)). The four lengths share the design but not the
  communities, so the columns are not paired.
- **Composition** (test samples without a host, median errors): pe and se unchanged (species' depth −3.8% and −2.6%).
  ONT −6.1% (v27 −5.5%). PacBio still overestimates its species' depth, at +27.6% (v27 +21.0%), with explained share
  +0.17 and unknown share −0.18. This has been open since v19.
- **Ancestry sites of the errors' reads:**
  - The AUCs of FN own reads against FP genus reads are lower than v27's: pe 0.575 → 0.556, se 0.583 → 0.546, pb
    0.467 → 0.395, ont 0.510 → 0.493.
  - The ancestry index changed to version 2 (`3ef6ad6`) and the error mix changed with the models, so this report does
    not separate the two.
  - The ancestry features carry no measurable F1
    ([ancestry true-positive test](../2026-10-09-ancestry-true-positive-test/README.md)), so this does not touch the
    calls.
- **The build audit's fixes show:**
  - `c_gates*.tsv` count the unseen species (G2/G3), with a table per model kind (G9).
  - Parity checks the gate models and both kinds (G21): all equal.
  - The marker-truth F1 counts no held-out species as missed.

## Recommendations

1. **Make GBM the database's default kind again** (`build_gtdb_database.py --model gbm,forest`, the trainer's default
   `--model gbm`). It is better in every read type at r226, by 0.005-0.009 on the scenarios. Keep the forest in the
   database for `--model_kind forest`. Re-check the real mocks with GBM (D6331's pe misses, and the gut row above),
   because that is where the forest looked better.
2. **Drop `--model-gates` from the default, or reduce it to one extra gate per read type.** The models are
   interchangeable across gates, and they cost ~70 minutes of the 4:16 build. What they did give is the deep samples
   at other gates (section 3). If that check is wanted, profile only the test collection at those gates: those rows
   are 28% of the extra cost.
   - A cheaper check of what the gate models add: score the C 1 test rows with the database gate's model (the
     trainer's `--extra-test-file`, as v27 did) in one build.
3. **Merge `f1-integration`'s knob change** (E2) before the next build, so that a 0.0001 margin cannot reset a model
   to 0.5 (7 of 24 here). Which samples the knob should serve, the design's or the held-out species', is the open
   choice behind it.
4. **Keep the database gates `pe:4,se:4,pb:16,ont:8` by stride for now.**
   - se C 8 and ont C 16 are free on the test rows and on the deep scenario samples, and save another 17% and 11% of
     the alignment CPU. They are the next step, once se C 8's forest model has its knob (recommendation 3).
   - The test-only rows at C 32 and at pe/se C 16 can go: they lose F1, as v27 already suggested.
5. **Label the grid rows' CPU in the console as the emulation's**, or leave them out. A grid sparse index's speed and
   memory need a real sparse-index run (`measure_performance.sh --test_sparse`).
6. Open, not started here: ONT's rowless species (3.7%), PacBio's depth bias (+27.6%).
