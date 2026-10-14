# Sparse gates on the real mock samples: ungated, an emulated hash sparse index, and the grid's sparse index at C 4

2026-10-14, protal 0.8.4 at `5369c99` (the latest commit when started), WSL on 4 cores (`taskset -c 0-3`; the machine
shared with other jobs, load 2-6).

**The question** (the user's): on the local mock samples, how do these compare?
- **`C0`:** an ungated run.
- **`C4s`:** `--sparse_first 4 --sparse_seeding`. The reads kept are seeded from the hash pick's cores alone, as an
  index of only those would seed them; the whole index is loaded.
- **`C4g`:** `--sparse_first 4 --sparse_first_pick grid`. The grid's sparse index, of grid level 2 and up, is loaded:
  4.13 M of the index's 16.38 M entries, on 3.70 M of its 10.25 M cores.

**Also run:** `dbgates`, the database's own gates (pe 4 stride, HiFi 16, ONT 8) on the whole index, the run
`measure_performance.sh` times.

**Data.**
- **Database:** the reduced r226 database of the [mock local test](../2026-10-13-mock-local-test/README.md) (761
  species), its inputs made again and built by `5369c99` with the local test's settings (pe 150 bp, HiFi, ONT; seed 1).
  - **Models at more gates:** C 1 and C 4 for every read type, besides the database's gates (`--model-gates
    pe:1/4,pb:1/4,ont:1/4`). So `C0` is scored by models trained ungated, and `C4s`/`C4g` by models trained at C 4.
  - **Model kinds:** random forests, protal's default since 0.8.4, and gradient-boosted trees.
  - **Size:** built in 23 min; `database.protal` 158 MB.
  - **Its own test sets:** F1 0.960-0.967 at every gate.
- **Samples:** the twelve mock samples in `~/Mocks` (1M pairs or ~1 Gb each).

**Method.**
- **Timing:** [measure.sh](measure.sh) runs `measure_performance.sh --test_sparse` per sample with the levels `0 4s 4g`,
  one run each, on 4 threads, without strain MSAs.
- **Scoring:** [collect.py](collect.py) copies each level's profiles into a run folder of the mock benchmark and scores
  them against the truth tables (`score_mocks.py score`).
- **Order:** [run_all.sh](run_all.sh) does everything in order.
- **GBM runs:** [gbm_runs.sh](gbm_runs.sh) does the same with the GBM models at `dbgates`, `C0` and `C4g`
  (`score_mocks.py run --protal-args "--model_kind gbm ..."`).

## Results

### Calls and abundances

Over the 10 samples with a truth of species, the host control left out, with the default forest models
([results/scores.tsv](results/scores.tsv), [results/samples_by_run_all.tsv](results/samples_by_run_all.tsv)):

| Level | TP | FP | FN | F1 | pe TP/FP/FN | long reads TP/FP/FN | Mean L1 (7 abundance samples) |
|---|---|---|---|---|---|---|---|
| `dbgates` | 197 | 9 | 80 | 0.816 | 149/0/76 | 48/9/4 | 0.512 |
| `C0` | 198 | 9 | 79 | 0.818 | 150/1/75 | 48/8/4 | 0.511 |
| `C4s` | 202 | 7 | 75 | 0.831 | 154/0/71 | 48/7/4 | 0.512 |
| `C4g` | 200 | 8 | 77 | 0.825 | 152/0/73 | 48/8/4 | 0.506 |

The four levels call within a few species of each other. The differences per sample:
- **D6331 DNA (pe):** both sparse levels call 11 of 15, the others 10.
- **hCom2 (pe, closed):** members called: `C4s` 56, `C0` and `C4g` 54, `dbgates` 53. No false positive in any.
- **Treichel (pe):** both sparse levels 63/0/0; `C0` 62/1/1, `dbgates` 62/0/1.
- **The long reads' false positives** are *Escherichia* and *Salmonella* congeners at low abundance (see the model
  kinds below):
  - ATCC HiFi: 4 at `dbgates`, 3 at `C0`, 1 at `C4s`, 2 at `C4g`;
  - D6331 ONT: 2 at `dbgates` and `C0`, 4 at `C4s`, 3 at `C4g`.
- **Unchanged:** the negative controls (no calls at any level), D6331 HKU and Auckland (12 of 15 at every level), and
  the strain mix.

**Abundances.**
- **`C4g` stays closest:** Bray-Curtis 0.0002-0.025 to the `dbgates` profile; in the sample where the levels differ
  most, D6331 DNA, `C0` is at 0.037 and `C4s` at 0.053.
- **L1 to the truth:** `C4g`'s is lower than ungated in 5 of the 7 abundance samples (HKU 0.185 against 0.189, D6331
  DNA 0.918 against 0.940) and within 0.002 in the other two.
- **`C4s` shifts them most:** Bray-Curtis up to 0.053 on paired-end reads, and HKU's L1 rises to 0.207.

### Records, time and memory

Per sample and level in [results/resources.tsv](results/resources.tsv).

| | `dbgates` | `C0` | `C4s` | `C4g` |
|---|---|---|---|---|
| pe records against `C0` | −0.1% | | −4.4% to −17.1% | −0.7% to +1.5% |
| Long-read records against `C0` | −0.2% | | −0.0% to −0.7% | −0.0% to −0.4% |
| Wall time, all 12 samples | 67 s | 116 s | 85 s | 72 s |
| Wall time, the 4 long-read samples | 41 s | 86 s | 47 s | 44 s |
| Peak memory | 0.6 GB | 0.6 GB | 3.4-3.5 GB | 0.5 GB |

- **Records:** the grid pick seeds from 1/4 of each gene's cores and anchors a short read's gene on one seed, and keeps
  the ungated run's records. The hash pick's sparse seeding loses 4-17% of the paired-end records.
- **Time:** ungated long reads take twice as long, as at r226. `C4g` runs as fast as the database's gates.
  - **Paired-end samples:** they take 2-7 s each, about a second of it the index load, so their differences are noise
    (one run each).
- **Memory:** this database is too small to show a sparse index's saving.
  - **The index:** 0.1 GB of values; 0.2 GB key map, the compact one, protal's default since `ab6bf25`.
  - **`C4g`:** its peak, 0.5 GB, is 0.1 GB below the others'.
  - **`C4s`'s 3.4 GB:** `--sparse_seeding` loads the index with the direct key map (3.2 GB), to report what an index of
    1/C of the cores would hold, and only then makes the compact one. Its memory and load time (1.9-4.7 s against
    0.5-1.0 s) are the emulation's, not those of a sparse index.

### A leak of the grid gate on a sample with host DNA

Treichel's mock is in mouse DNA (its unknown share 0.76).

| | `C4s` (hash) | `dbgates` (stride) | `C4g` (grid) |
|---|---|---|---|
| Pairs given up by the gate | 97% | 77% | 49% |

**What the grid gate costs there** (`logs/pe1.sparse4g.log`):
- **The gate:** it looks up every core of a read and counts the grid's entries, 4.4 sampled lookups a pair against 2.0
  for the hash. A foreign read meets a grid entry more often.
- **The anchors:** a short read's gene anchors on a single grid seed. 992 k mates got an anchor, against 599 k ungated.
- **The candidates:** 2.48 M candidate alignments were tried, against 0.75 M ungated; the k-mer screen refused 2.03 M of
  them.
- **The time:** the alignment took 3.0 s against 2.5 s ungated and 1.7 s for `C4s`.
- **The output:** the same records (15,266 against 15,044) and the same calls.

In the D6331 samples without a host the grid gate gives up as much as the hash gate (93-95%). Host-rich samples (human
gut, mouse) at r226 should be timed with the grid pick before it is recommended for them.

### The model kinds

The gate comparison holds with the GBM models. The model kind, though, moves the long reads' calls more than any gate
does:

| Level and model | TP | FP | FN | F1 | pe TP/FP/FN | long reads TP/FP/FN | Mean L1 |
|---|---|---|---|---|---|---|---|
| `dbgates`, forest | 197 | 9 | 80 | 0.816 | 149/0/76 | 48/9/4 | 0.512 |
| `dbgates`, GBM | 192 | 2 | 85 | 0.815 | 144/1/81 | 48/1/4 | 0.567 |
| `C0`, GBM | 193 | 2 | 84 | 0.818 | 146/1/79 | 47/1/5 | 0.568 |
| `C4g`, GBM | 194 | 1 | 83 | 0.822 | 146/0/79 | 48/1/4 | 0.567 |

- **The forest (0.8.4's default):** it calls the long reads' *Escherichia* congeners. With it, D6331 HiFi is 14/3/1 and
  D6331 ONT 14/2/1; with GBM 14/0/1 for both, and ATCC HiFi 19/4/1 against 19/1/1.
- **GBM:**
  - **D6331 Auckland (pe):** it misses *Veillonella rogosae* (20% of the cells, 15× deep, P 0.34) and *E. coli* (8.7%),
    which the forest calls; its L1 rises from 0.41 to 0.85.
  - **hCom2:** it calls 4 fewer of the members.
- **Overall:** F1 is the same for both kinds (0.815-0.825). Neither kind is the better one on these samples.

## What this means

- **Calls:** at C 4 on these real samples, the grid's sparse index calls as the ungated run and the database's stride
  gate do, and keeps the most accurate abundances. The hash pick's sparse seeding loses 4-17% of the paired-end records
  and shifts the abundances most.
- **Speed:** the grid runs at the database gates' speed, twice as fast as ungated on long reads.
- **Not measurable here:** the memory a sparse index saves needs a database of r226's size. To run on the cluster:
  `measure_performance.sh --test_sparse` with `SPARSE_SHORT="0 4 4g"` on r226.
- **Not measured at all:** a host-rich sample's cost with the grid pick at r226.
- **For the model default:** on these samples the forest's long-read false positives and the GBM's missed abundant
  paired-end species weigh about the same; the mock benchmark at r226 should decide.

## Files

| File | What it does |
|---|---|
| [run_all.sh](run_all.sh) | the clone, the inputs, the database, the measurements and the scores, in order |
| [measure.sh](measure.sh) | `measure_performance.sh --test_sparse` (levels `0 4s 4g`) per mock sample |
| [collect.py](collect.py) | the levels' resources side by side, their profiles into run folders, the scores |
| [gbm_runs.sh](gbm_runs.sh) | the GBM models at `dbgates`, `C0` and `C4g`, and the scores of all seven runs |
| [results/resources.tsv](results/resources.tsv) | per sample and level: time, memory, records, reads given up, species against `dbgates` |
| [results/scores.tsv](results/scores.tsv) | per sample and level (forest): calls, F1, L1 and Bray-Curtis against the truth |
| [results/samples_by_run_all.tsv](results/samples_by_run_all.tsv), [results/summary_by_run_all.tsv](results/summary_by_run_all.tsv) | all seven runs' scores (`score_mocks.py`) |

**Kept in WSL:**
- `~/c4test`: the clone (built), the inputs, the database and the measurement folders, with their logs.
- `~/Mocks/0.8.4_5369c99_*`: the runs' profiles.
