# The r226 v24 build against v23: the gain is the easier mix, and what the indexed alleles did

**Data.**
- **v24:** SLURM 24161358, q512n11, 2026-10-11 17:42 to 19:50, 84 threads. protal and the simulator at `7d5df4b`
  (0.8.1), the scripts at `7d5df4b` at the start and `2caffd7` at the end (`build_metadata.tsv`), seed 1,
  `--index-alleles 0.01` given on the command line. The share archive `local/protal0.8.1_r226_v24_share.tar.gz`
  (543 MB) is unpacked to `local/v24`.
- **v23** ([2026-10-11-r226-v23](../2026-10-11-r226-v23/README.md), `local/v23`, `8d8c89b`): the same 47,403
  genomes of 34,956 species, the same hold-out (12,142 species), design, scenarios and seed, and the same
  `--index-alleles 0.01`.
- **Commands:** [builds_by_class.py](builds_by_class.py) (v23's script with the build directories found by a glob),
  run in WSL on cores 0-3 as `builds_by_class.py local <read type> v23 v24` for pe, se, pb and ont;
  output [builds_by_class.txt](builds_by_class.txt). The index and allele numbers come from
  `logs/training_db_index.log` and `logs/index_and_package.log` of both builds.

**Questions (the user's, 2026-10-11).** Is v24 really better than v23? `--index-alleles 0.01` was on by accident:
is there any sign that it mattered?

## Summary

- **v24's headline gain is the easier sample mix, not better error rates.** With species held out, v24 is ahead by
  pe +0.0013, se +0.0011, pb +0.0003, ont +0.0001. Applied to the same class mix, v24's error rates are worth
  pe +0.0003 / −0.0005, se +0.0001 / −0.0009, pb −0.0003 / −0.0016 and ont −0.0011 / −0.0022 (at v24's mix / at
  v23's). The independent test sets move both ways within noise (pe −0.0012, se +0.0015, pb +0.0050, ont −0.0029).
- **Why the mix is easier:** `--insilico-strains auto` cut the in-silico strains from 24,727 to 8,575 of the
  one-genome species, and from 5,136 to 2,497 of the species with real strains. No new genomes were downloaded, so
  the real strains stayed at 14.5% of the simulated species. Strains fell from 53.3% to 29.0% of the simulated
  species, and representatives (the exact database genome, the easiest class) rose from 31% to 62% of the present
  taxa.
- **The rates did move, in opposite directions.** False positives beside a held-out congener fell by 8-10% in every
  read type. In-silico strains are missed 22-39% more often, and real strains without alleles 6-23% more often. A
  model trained on fewer strains reads divergent evidence more often as a missing congener's reads. This is a
  trade-off between the two error classes, not a gain on both.
- **Real samples are mostly strains, not representatives,** so v23's mix is the closer of the two to real use. On
  that mix v24's model is slightly worse.
- **The indexed alleles cannot explain any v23 → v24 difference.** `--index-alleles 0.01` was on in v22, v23 and
  v24. They also cannot be measured from v24: every run of the build used them. The one isolated measurement, the
  strain false-negative study, found no effect: `a16` → `a16i0` gained 4 strains and lost 2 of 600 (p = 0.69).
- **The index grew and costs more now.** The training database has 118.0 M allele seeds (v23 101.6 M; 16 alleles).
  The finished database has 202.5 M (v23 114.8 M), because its alleles now come from every genome (`f4850d1`). The
  models were therefore trained with 1.7× fewer indexed allele seeds than the database they ship in. Scaled from
  v23's measured cost, the seeds take ~1.8 GB of `database.protal` and ~1.2 GB of profiling memory.
- **The weighted site shift was off in v24** (protal `7d5df4b`; `353b88b` made it the default later). Profile v24
  with `--no_weighted_site_shift` on a protal at `353b88b` or later. protal's help says this applies to databases
  "built before 2026-10-11", but v24 was built on that day and still needs it. The r226 test of the shift is still
  to come.

## 1. The headline numbers

F1 at the knob protal calls at, from each build's `console.log` (test-set knobs: pe 0.74 in both builds; ont 0.53
in v24, 0.5 in v23; the others 0.5):

| | pe held out | pe test | se held out | se test | pb held out | pb test | ont held out | ont test |
|---|---|---|---|---|---|---|---|---|
| v23 | 0.9595 | 0.9709 | 0.9650 | 0.9662 | 0.9791 | 0.9828 | 0.9741 | 0.9820 |
| v24 | 0.9608 | 0.9697 | 0.9661 | 0.9677 | 0.9794 | 0.9878 | 0.9742 | 0.9791 |
| change | +0.0013 | −0.0012 | +0.0011 | +0.0015 | +0.0003 | +0.0050 | +0.0001 | −0.0029 |

Species held out, pe: FP 4,049 → 3,780 and FN 2,255 → 2,281. The gain is all in the false positives. The test sets
are small (pb 2,451 taxa with 14 FP and 21 FN), and v23's bootstrap intervals on such sets were ±0.002-0.005.

The scenarios (hold-in, species held out) gained up to 0.003 on the complex communities (pe moderate 0.972 →
0.975, soil 0.956 → 0.958, soil_shallow 0.948 → 0.951). Their FP are the class that fell. Gut lost (pe 0.991 →
0.989, ont 0.994 → 0.988). The 4-sample hold-out sets are unchanged within ±0.002.

Composition, as the median over the test samples without a host: se's unknown share −0.056 → −0.025, ont's species
depth −8.1% → +2.8%, PacBio's explained share +0.196 → +0.183. These depend on the mix too.

## 2. Mix against rates

F1 at knob 0.5 with species held out, computed from the class error rates of one build (row) and the class shares
of the other (column). The script uses the classes of the v23 report:

| F1 from the class rates | pe v23 mix | pe v24 mix | se v23 | se v24 | pb v23 | pb v24 | ont v23 | ont v24 |
|---|---|---|---|---|---|---|---|---|
| v23 rates | 0.9595 | 0.9605 | 0.9650 | 0.9660 | 0.9791 | 0.9797 | 0.9741 | 0.9753 |
| v24 rates | 0.9590 | 0.9608 | 0.9641 | 0.9661 | 0.9775 | 0.9794 | 0.9719 | 0.9742 |

The mix alone is worth +0.0006 to +0.0012 at v23's rates and +0.0018 to +0.0023 at v24's. v24's rates add
between +0.0003 and −0.0022.

Present taxa (pe, species held out): representatives 31.1% → 62.4%, in-silico strains 53.8% → 19.4%, real strains
15.1% → 18.1%. The share of real strains with alleles among real strains is unchanged (57.0% → 56.9%), so the
allele classes are comparable between the builds.

Error rates by class, species held out, v23 → v24:

| class | pe | se | pb | ont |
|---|---|---|---|---|
| FP beside a held-out congener | 1.41% → 1.29% | 1.32% → 1.21% | 1.79% → 1.61% | 1.82% → 1.66% |
| FP beside present congeners only | 0.61% → 0.59% | 0.88% → 0.85% | 1.65% → 1.54% | 1.21% → 1.35% |
| FP no congener in the sample | 1.22% → 1.24% | 1.81% → 1.85% | 11.2% → 11.8% | 5.88% → 7.32% |
| FN in-silico strain | 2.75% → 3.35% | 2.79% → 3.40% | 1.88% → 2.62% | 2.54% → 3.35% |
| FN real strain, no alleles | 4.88% → 6.01% | 5.43% → 6.28% | 7.25% → 7.97% | 7.06% → 7.46% |
| FN real strain, alleles | 6.51% → 6.51% | 6.40% → 6.38% | 6.74% → 7.04% | 7.31% → 7.42% |
| FN representative | 1.84% → 1.90% | 1.78% → 1.81% | 1.07% → 0.99% | 1.51% → 1.64% |

The two big moves are each several binomial standard errors (pe in-silico strains 1,137 of 41,407 against 499 of
14,888; FP beside a held-out congener 3,104 of 220,460 against 2,799 of 217,657). Correlation within samples makes
that an overstatement, but the direction is the same in all four read types.

The real strains with alleles did not gain from 16 alleles instead of 4 (pe 6.51% in both). The strain
false-negative study expected a few points there. The rise of the other strain classes, which is plausibly an
effect of training on fewer strains, may hide such a gain. A refit of v24's tables with the strain rows weighted
up to v23's share would separate the two; it was not run.

## 3. The indexed alleles

| | v23 training DB | v24 training DB | v23 finished DB | v24 finished DB |
|---|---|---|---|---|
| strain alleles (per species and gene) | 3.70 M (4) | 5.03 M (16) | 4.22 M (4, half the genomes) | 9.23 M (16, all genomes) |
| alleles indexed (>= 0.01 from the representative) | 1.33 M | 1.93 M | 1.50 M | 3.56 M |
| allele seeds in the index | 101.6 M | 118.0 M | 114.8 M | 202.5 M |
| `database.protal` | 24.7 GB | 25.5 GB | 25.3 GB | 26.7 GB |

- **On in every compared build.** v22, v23 and v24 all ran `--index-alleles 0.01`, so no build-to-build
  difference can come from it. Within v24 every training and test run used the indexed training database, so
  nothing in this build measures it.
- **The isolated measurements say it does not matter for F1:** the strain false-negative study
  ([2026-10-11-strain-false-negatives](../2026-10-11-strain-false-negatives/README.md)) gained 4 and lost 2 strains
  of 600 with 16 indexed alleles (p = 0.69), and aligned pairs did not change. The ancestry truth test measured
  +0.003 for the default model.
- **Training does not match the finished database.** The models saw 118 M allele seeds and are used with 202.5 M.
  The unindexed alleles have the same gap by design (`f4850d1`). The index makes the gap larger and has not been
  tested at r226.
- **Cost, scaled from v23's 114.8 M seeds** (database +1.0 GB for the seeds alone, profiling index ~0.7 GB, build
  memory after pass 2 +2.8 GB): for 202.5 M, ~1.8 GB of `database.protal`, ~1.2 GB of profiling memory and ~4.9 GB
  of build memory. These are estimates, not measurements. `protal_db`'s build peak was 46.3 GB, against 44.9 GB in
  v23.

## 4. What to do

1. **Do not count v24's +0.001 as a model gain.** Compare builds at one class mix (`builds_by_class.py`), or keep
   the mix fixed between builds.
2. **Keep strains near half of the simulated species.** Run the larger download (`download_gtdb.py --species 40000
   --per_species 4`) before `--insilico-strains auto` lowers the in-silico share, or give `--insilico-strains` a
   share that keeps the total strain share near v23's.
3. **Build the next database without `--index-alleles`** (the default is off). This saves ~1.8 GB of disk and
   ~1.2 GB of profiling memory, removes the training/finished index gap, and by every isolated measurement costs
   no F1.
4. **Profile v24 with `--no_weighted_site_shift`.** Its models were trained without the shift. The next build
   tests the shift at r226 (it is on by default in the builds since `353b88b`).
