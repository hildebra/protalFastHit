# A sparse index by the genes' grid: proposal and implementation

Date: 2026-10-13. Branch `audit-fixes`, on `0ab85f4` plus the uncommitted stride default
([r226 v27 gates](2026-10-13-r226-v27-c-gates.md), follow-up section). The proposal (no code, no run: a design, with
estimates from earlier reports), then its implementation and first measurements on the mini database (follow-up).

## The question

A gated run's reads are looked up first by 1/C of their k-mer cores. The user asked whether the index could then hold
only the keys those lookups need, to save memory. The stride pick (every C-th core of a read in read order, the
default since 2026-10-13) is chosen on the read, so any key may be one it picks: the whole index must be loaded. The
hash pick can be served by an index of 1/C of the keys ([sparse seeding](2026-10-13-sparse-seeding/README.md)), but
the hash lost 2.5-17 times more present species than the stride at r226.

The user's idea: choose the sampled cores at build time by their place in the original gene.

## The idea: a grid on each gene

The build keeps, of each gene's cores in gene order (its syncmers), every C-th: core ordinal i is kept at C if
i mod C = 0. The grids nest: a core kept at C = 8 is kept at 4 and 2. One file can serve every C if each occurrence
stores its level, ctz(i) capped at 5, and a run gated at C loads only the occurrences of level log2(C) or more.

The read cannot tell which of its cores lie on the grid, because its offset in the gene is unknown. So it looks up
every core, and only the grid ones hit:

- **The hits are evenly spaced.** A read covering a stretch of a gene hits that stretch's grid cores, exactly every
  C-th of its cores. That is the stride's spacing, with no hash gaps. A sequencing error loses only the cores it
  touches; it does not shift the grid, which is the reference's.
- **A miss costs one key-map probe.** At r226 a probe is ~20 of the 175 ns of a lookup; the flex cells take ~60%
  ([seeding profile](2026-10-12-seeding-profile/README.md)). A hit's block holds only the grid occurrences, about 1/C
  of today's, so it scans fewer flex cells.
- **The kept reads are seeded from their grid hits alone.** There is no second round of lookups. This is what
  `--sparse_seeding` emulates for the hash.

## What it would save at r226 (estimates)

The values are ~1/C of today's 26.6 GB, as for the hash: ~6.7 GB at C 4, ~3.3 GB at C 8. The keys do not shrink as
much. A core occurs ~30 times on average (2.9 G entries on ~0.1 G keys), and a key is kept if any of its occurrences
is on a grid: about 1 − (1 − 1/C)^30, ~98% of the keys at C 8. Today's key map is 3.2 GB whatever it holds, so that
matters only for a compact key map, which would save little here.

A HiFi run at C 8 would then need ~13 GB instead of ~36 GB, as with the hash's sparse index and today's key map.

## The catches

- **Strains shift the grid.** An indel or a SNP that adds or removes a syncmer early in a gene shifts the ordinals
  after it. Two strains' copies of a gene then have their grids at different cores. Each copy keeps its own grid,
  so a read of either strain still hits its own copy every C cores, and sometimes the other's too.
- **Conservation could choose the phase.** Instead of starting at the first core, the phase of each gene's grid
  could be chosen so that its kept cores are the most conserved in the species. That helps the gate on novel strains
  and makes copies' grids agree, but it biases the seeds away from variable regions, which carry the strain signal.
  It is a second step at most.
- **The level must be stored or recomputed.** An entry holds its gene and its position, not its ordinal. Either the
  build writes 3 level bits per entry (values +~6%, or ordered by level within a key's block), or the load recomputes
  every gene's syncmer ordinals from the packed genes (one pass over 14.5 M genes at r226).
- **The models must be trained with it.** As with sparse seeding, the seed-count features shrink by C. The build would
  train the models of a grid database on grid-seeded runs and record that in `build_settings.tsv`; with
  `--model-gates`, one model per gate.
- **Short reads see few grid cores.** A 150-base mate with ~20 cores has 2-3 grid cores at C 8. It is the same count
  as the stride's, and r226 v27 showed that pe/se hold their F1 up to C 4 with the full index. Whether seeding from
  those cores alone holds is what has to be measured.

## What to measure first

1. **Emulate it** as `--sparse_seeding` emulates the hash's index. Add a developer option `--grid_seeding C`: with
   the full index loaded, a lookup keeps only the occurrences whose ordinal in their gene is a multiple of C (the
   ordinals computed once at the start from the preloaded genes). The gate then counts grid hits, and the kept reads
   are seeded from them. The log says how many entries and keys such an index would hold per C.
2. **Compare locally** (the archaea and small real worlds), at the same C, for pe, HiFi and ONT:
   - the stride with the full index (today's default);
   - the hash with sparse seeding;
   - the grid.

   Measure the records against the ungated run, the reads given up, and the aligning time.
3. **At r226**, `measure_performance.sh --test_sparse` with grid levels, then a build whose models are trained grid
   seeded, at the gates it would ship.
4. **If the grid holds**, build it into the database: the level bits, a load of one level, and runs warned when gated
   otherwise.

The stride with the full index stays the default until the grid has been measured.

## Follow-up: implemented, first measurements (2026-10-13, evening)

The user asked to build it into the database directly, choosing the core each window keeps by the species' strain
alleles and the genes' conservation rather than by its ordinal. Code on `4dfc549` plus the uncommitted grid changes
(`src/SequenceUtils/GridLevels.h`, `Hash/Seedmap.h`, `Hash/KmerLookup.h`, `Core/ChainAnchorFinder.h`, `Build.h`,
`RunProtal.h`, `Options.h`, `Utilities/BuildSettings.h`; tests `tests/test_GridLevels.cpp` and an e2e test; scripts
`collect_training_data.py --grid_levels`, `build_gtdb_database.py --test-grid` and `--sparse-first-pick grid`,
`measure_performance.sh --test_sparse` levels such as `8g`). Measured on the mini database only.

### What was built

- **Levels by a tournament, not by ordinal.** Each gene copy's cores (its syncmers, in gene order) are split into
  windows of 2^l ordinals; each window keeps the better of its halves' winners, so a core's level is the last round it
  won (0-5) and the cores kept at C are a subset of those kept at C/2. A gene shorter than a window keeps its winner up
  to level 5, so every gene keeps a core at every gate (at C 32 the mini database keeps 1/29.8 of the cores, not 1/32).
- **The score of a core** (`GridLevels.h`): the share of the species' strain alleles (`strain_alleles.tsv`) covering
  its 31-mer that keep it a gate hit (no edit in the 15-base core, no indel in the 31-mer, at most two substitutions in
  the 16 flank bases, the gate's score 14), pooled with a prior from the conservation of its columns within the
  family's genera (`column_weights.tsv`, the `within` code, counted as two alleles; 0.5 without weights); minus
  0.5 × log2(entries of its key) / 11, so that a core many genes share loses; minus 0.05 per window length of distance
  from the window's middle, so that ties stay spread. A core whose key the index leaves out never wins. This answers
  the catch above: the grid is chosen by conservation, not phase, so strains' copies agree where their alleles do.
- **Storage.** The level takes three bits of the entry's position (bits 17-19; index feature `grid levels`); a
  database whose positions need more is built without levels. In memory the packed layout keeps the three bits only
  for a run that needs them (`Index in memory: ... grid level 3`), so other runs' index is as large as before.
- **The pick and the sparse index.** `--sparse_first_pick grid` looks up every core and counts the entries of level
  log2 C or more; the gate's hits and a kept read's seeds come from those alone. `--sparse_index auto` (default) loads
  only those entries when every sample of the run is grid-picked; `off` loads the whole index and filters per lookup.
  Both seed identically (the unit test and the e2e test compare records and reads given up; on the simulated sample
  below the SAMs agree to the byte apart from the `@PG` line).
- **The sparse load** decodes every chunk of the file, keeps the entries of the level, and writes the key map's
  pointers itself (a prefix sum per chunk: every key with entries in the file passed the build's limits on a key and a
  block, so a subset of its entries does too). A first version with a temporary key map per chunk, a walk of the keys
  one by one and `BuildValuePointers` took 2.55 s on the mini database against 0.11 s for its whole index, because the
  key map has 2^30 keys whatever the reference; now 0.25-0.44 s on 4 threads.

The mini database's build log:

```
Grid levels: 71054 cores, 99.858% of them covered by strain alleles, the others scored by the column weights; mean
conservation 0.589; kept at C 2: 35618 (1/2.0, conservation 0.610), C 4: 17902 (1/4.0, conservation 0.642), C 8: 9032
(1/7.9, conservation 0.656), C 16: 4600 (1/15.4, conservation 0.658), C 32: 2383 (1/29.8, conservation 0.658)
```

The kept cores' estimated survival in the species' strains rises from 0.589 for all cores to 0.642 at C 4 and 0.658
at C 8 and above: the choice does what it was meant to.

### Records on the mini database

Data: `scripts/mini_db/build_mini_db.sh` (its simulated release), two paired-end samples
(`simulate_metagenomes --genome_table genomes.tsv -n 2 --total_read_pairs 200000 --species_per_sample 20 --seed 7`)
and two HiFi samples of the same communities (`--long_samples`, 60 Mb each, `hifi:15000:3000:3`), each run once on 4
threads by `measure_performance.sh -t 4 --test_sparse` with `SPARSE_SHORT="1 2 2g 4 4g 4s 8 8g 8s 16 16g"` and
`SPARSE_LONG="1 4 4g 8 8g 8s 16 16g 16s 32 32g"` (HiFi with `PROTAL_ARGS=--no_profile`: the mini database has no
PacBio model). Plain levels are the stride on the whole index; `g` the grid on its sparse index; `s` the hash pick
with `--sparse_seeding`, an emulated index of the hash pick's 1/C of the cores. Records written against the ungated
run, the two samples:

| C | pe: stride, whole | pe: grid, sparse | pe: hash, sparse | HiFi: stride, whole | HiFi: grid, sparse | HiFi: hash, sparse |
|---|---|---|---|---|---|---|
| 2 | −0.02/−0.03% | −1.0/−1.1% | | | | |
| 4 | −0.15/−0.27% | −3.0/−3.5% | −5.8/−7.1% | 0/0 | −0.01/−0.03% | |
| 8 | −1.2/−1.8% | −8.1/−10.1% | −20.1/−22.4% | 0/0 | −0.04/−0.12% | −0.15/−0.26% |
| 16 | −3.0/−4.3% | −47/−48% | | 0/0 | −0.09/−1.0% | −0.22/−0.36% |
| 32 | | | | 0/0 | −0.34/−2.0% | |

- **The gate itself** gives up as many pairs with the grid as with the stride (pe C 4: 42.6% against 42.5%; C 8:
  43.2% both), and fewer than the hash (43.7%, 46.4%): the grid's hits are as evenly spaced as the stride's.
- **The records lost are the sparse index's.** A kept pair seeded from its grid cores alone has a quarter of its
  seeds at C 4. Paired-end reads lose 3% of their records at C 4, half of what the hash pick's 1/C index loses, and
  are unusable at C 16 (about one grid core per mate). The unknown share of the profiles grows with them (pe1:
  0.0067 ungated, 0.0136 at 4g, 0.0288 at 4s); the mini world's three called species did not change.
- **HiFi reads** have hundreds of cores each and lose at most 0.1% up to C 8; at C 16 and 32 one sample lost 1-2%,
  more than the hash's index there.
- **Times and memory** say nothing here: the runs take 1-3 s and the values of the mini index are kilobytes beside
  the 3.2 GB key map (every run peaked at 3.1-3.2 GB).

### What is open

1. **r226.** A build with `--test-c-gates` now also profiles the test set with the grid pick at each gate with a model
   of 2 or more (`--test-grid auto`; rows `<sample>__c<C>g` in `c_gates.tsv`), scored by the stride-trained models of
   the nearest gate: the first F1 of the grid. Then `measure_performance.sh --test_sparse` with `4g 8g 16g` for the
   memory and load time of the sparse index at r226 (estimated ~6.7 GB of values at C 4, ~3.3 GB at C 8).
2. **Models trained grid-picked.** If the F1 holds, `build_gtdb_database.py --sparse-first-pick grid` builds databases
   whose models are trained on grid-picked runs (fewer seeds change the seed-count features).
3. **Paired-end seeding from a denser level than the gate.** Gating at C 4 by level 2 while seeding kept reads from
   level 1 (an index of 1/2) would trade memory for paired-end records (the 2g row: −1% instead of −3%). Not built.
4. **The sparse load still decodes the whole file.** Its time is about the whole load's decode plus a filter; at r226
   it saves memory, not load time.

## Follow-up 2: the best core per window, or a stride on the gene? (C 4)

The user asked how the grid differs from the stride, and whether the stride's rule could run on the genes. It can: a
stride on the gene (every C-th core, C apart) is what the grid's tournament does when every score is equal (ties go
to the window's middle: of 8 cores, C 2 keeps 0, 2, 4, 6 and C 4 keeps 2 and 6). The scored grid lets each window keep
its best core instead, anywhere in the window, so its kept cores are 1 to 2C − 1 apart. Which costs fewer records was
measured at C 4 only, as the user asked.

**Variants** (`protal --build --grid_choice`, a developer option): `scored` (the default above); `even` (equal scores:
a stride on the gene); `phase` (every C-th core too, C apart, with the phase per gene and level chosen by the summed
scores of the cores it keeps, nested in the level below's: conservation chooses where the stride falls, not where
each core falls). The window's central half, first proposed as the middle variant, does not nest (a window's core must
be one of its halves'), so `phase` took its place.

**Data and commands.** The working tree of follow-up 1 plus `--grid_choice` (unit tests: 10 `Grid*` tests pass). Three
mini databases from one simulated release (`PROTAL_BUILD_ARGS="--grid_choice <c>" scripts/mini_db/build_mini_db.sh`;
genome table and full reference identical, 71054 cores each); 4 paired-end samples (`simulate_metagenomes -n 4
--total_read_pairs 200000 --species_per_sample 20 --seed 7`) and 2 HiFi samples of samples 1 and 2's communities (60 Mb
each, `hifi:15000:3000:3`); `measure_performance.sh -t 4 --test_sparse` with `SPARSE_SHORT=4g SPARSE_LONG=4g`
(`PROTAL_ARGS=--no_profile`, `REPEATS=1`), on the scored database also `4` and `4s`.

| | kept cores' conservation (C 4) | grid hits, pe1-pe4 (M) | pe records lost, pe1/pe2/pe3/pe4 | mean | HiFi records lost |
|---|---|---|---|---|---|
| stride, whole index (`4`) | | | 0.15/0.27/0.24/0.26% | 0.23% | |
| scored (`4g`) | 0.642 | 1.025/1.061/0.963/1.014 | 2.95/3.47/3.03/3.39% | 3.21% | 0.01/0.03% |
| even (`4g`) | 0.590 | 1.009/1.025/0.883/0.955 | 3.40/4.34/3.64/4.11% | 3.87% | 0.01/0.03% |
| phase (`4g`) | 0.598 | 1.021/1.047/0.896/0.976 | 3.68/4.51/3.98/4.30% | 4.12% | 0.10/0.17% |
| hash, emulated sparse index (`4s`) | | | 5.76/7.07/6.46/6.89% | 6.55% | |

- **The scored choice wins on every sample.** It keeps cores that the reads carry more often (0.642 against 0.590
  estimated survival in the species' strains), so a run finds 1.6-9% more grid hits, and loses 0.66 points fewer of its
  records than the gene stride (3.21% against 3.87%). The uneven spacing costs less than the better cores gain.
- **Choosing only the phase does not help.** One phase per gene and level cannot follow the conservation along the gene
  (0.598 kept), and it lost slightly more than the plain gene stride, also on HiFi (0.10-0.17% against 0.01-0.03%). Why
  it is below `even` was not looked into; one difference is that a phase keeps a core the index leaves out where the
  tournament gives the window to its other half's core.
- **The gate is unaffected**: all three give up 42.3-43.6% of the pairs, as the stride does (42.2-43.3%).
- **None of them closes the gap to the stride on the whole index** (0.23%): that gap is the sparse index's seeding
  from 1/C of a kept read's cores, not the choice of cores.

`scored` stays the default; `even` and `phase` stay as a developer option (`--grid_choice`) for a check at r226.

## Follow-up 3: single-seed anchors (C 4)

**Why the grid lost records.** A short read's mate makes an anchor only from two or more seeds of one gene on one
diagonal (`ChainAnchorFinder::FindPairs`, `FindAnchorsSingleRef`); a mate whose partner anchored on a gene may anchor
there on one seed (`RecoverAnchors`). Once anchored, the whole mate is aligned against the gene, so the seeds no longer
matter. A grid-seeded mate has about 1/C of its seeds (a 150-base mate ~23 cores, ~5-6 grid cores at C 4 if wholly in
the gene, fewer at a gene's end or where errors and strain SNPs hit grid cores), and the two-seed rule, a guard set
for reads with 20-odd seeds, drops the genes with one. The earlier runs show it growing with C: of the seeds, 99.6%
shared their gene with another at stride C 8 (whole index), 93.5% at grid C 8, 56% at grid C 16 (47% fewer records).

**The change.** With the grid pick, a short read's gene among its seeds that has no anchor anchors on its first seed
(`AnchorLoneSeeds`: `ExtractAndExtendAnchorFromSeed`, as the mate rescue does): a grid seed is an exact 15-base core
with its flanks matched on its gene's own grid, and fixes the read's place on the gene by itself. These anchors rank by
their extended length with the others, so only a mate without better anchors aligns on them (`--align_top`). Long
reads (chunks of hundreds of cores) are left as they were. `--no_single_seed_anchors` turns it off for comparison; the
seeding line counts them ("A anchors (S of them of a single grid seed)"). Also fixed: the mate rescue
(`KmerLookupSM::RecoverFromLookup`) read every entry of a grid-picked lookup's block on the whole index, where a sparse
index holds only the level's; it now skips the others, so `--sparse_index off` rescues as `auto` does.

**Data and commands.** The working tree of follow-up 2 plus this change; unit tests (33 of the seed, anchor, lookup
and grid suites) and the 11 sparse-first e2e tests pass (the e2e test checks the count and that
`--no_single_seed_anchors` writes fewer records). The scored mini database (`scripts/mini_db/build_mini_db.sh`), the 4
paired-end samples of follow-up 2 (`simulate_metagenomes -n 4 --total_read_pairs 200000 --species_per_sample 20 --seed
7`, the same reads), `measure_performance.sh -t 4 --test_sparse` with `SPARSE_SHORT="4 4g"`, and again with
`SPARSE_SHORT=4g SPARSE_ARGS=--no_single_seed_anchors` (`REPEATS=1`, profiles on). The candidate alignments and
anchors are the runs' log lines.

| pe1/pe2/pe3/pe4 against the ungated run | records written | mean | candidate alignments tried | profiles: Bray-Curtis, unknown share |
|---|---|---|---|---|
| stride, whole index (`4`) | −0.15/−0.27/−0.24/−0.26% | −0.23% | −0.1 to −0.2% | ≤ 0.0005; unchanged ±0.001 |
| grid, two-seed anchors (`4g --no_single_seed_anchors`) | −2.95/−3.47/−3.03/−3.39% | −3.21% | −9.9/−17.8/−6.5/−13.6% | 0.0028-0.0043; 0.0067 → 0.0136 (pe1) |
| grid, single-seed anchors (`4g`) | −0.04/−0.12/+0.03/−0.08% | −0.05% | +4.7/+4.7/+1.3/+3.1% | ≤ 0.0009; 0.0067 → 0.0072 (pe1) |

- **The records are back**, to within 0.12% of the ungated run on every sample, closer than the stride on the whole
  index. One sample writes 0.03% more than ungated: a gene with one seed in a full-index read makes no anchor there
  either, and now does.
- **The single-seed anchors** are 9.5-23% of the grid runs' anchors (20,376-65,812 per sample).
- **The cost** is 1-5% more candidate alignments than the ungated run (8-27% more than the two-seed grid run); the
  aligning times of these 1-3 s runs are noise. The lookups stay those of the grid.
- **Calls**: the mini world's three called species are the same in every run, so false positives from single-seed
  candidates to congeners cannot show here; that, the F1 and the cost at scale are for r226 (`--test-grid`).

## Follow-up 4: two-seed anchors on the whole index

The user asked that runs with the full index keep the two-seed anchors. Since single-seed anchors only pay where a read
has 1/C of its seeds, the grid now seeds by what is in memory:

- **Sparse index** (`--sparse_index auto` with every sample grid-picked): as in follow-up 3, a read kept is seeded from
  its grid entries alone, a short read's gene anchoring on one seed too.
- **Whole index** (`--sparse_index off`, or a run whose samples are not all grid-picked): the grid only gates. The
  gate keeps its lookups of every core, those without an entry of its level too (scored as misses), and a read kept is
  seeded from them with every entry counting (`ChainAnchorFinder::KeepSparseLookups`), so no core is looked up twice;
  anchors need two seeds of a gene, as in every other run.
- **`--sparse_index emulate`**: the whole index, seeded as the sparse index (the same records as `auto`). Every
  collection run of `collect_training_data.py` gets it where protal has it, so that the build's `--test-grid` rows,
  profiled in mixed runs on the whole index, and the training rows of a database picked by grid are seeded as a
  grid-picked run on a sparse index is. It changes nothing for samples of the other picks.

The run's log says which ("Grid pick: ..."). The e2e test now checks `auto` = `emulate` (records, pairs given up), and
that `off` gives up the same pairs, has no single-seed anchors and more than twice the seeds. Not yet measured: how the
grid gate on the whole index compares with the stride (its records should be close to the stride's at the same C,
the gate giving up as many pairs).

## Follow-up 5: held-out real strains (C 4)

The user asked whether the gates miss divergent strains. The test bed is the strain world of
[2026-10-11-strain-false-negatives](2026-10-11-strain-false-negatives/README.md) in WSL `~/strainfn`: GTDB r226 marker
genes of 30 focal species, 532 real test strains held out of the database (easy ≥ 99.5% ANI to the representative,
198; medium 98.5-99.5%, 248; hard < 98.5%, 86), 20 deep samples (60 pairs per strain) and 20 shallow ones (1-32 pairs
per strain, 600 instances), and the same presence model for every run (the small real world's paired-end GBM, seed 1,
trained on ungated runs with full seeding).

**Commands.** The tree of `19ccc0e` plus the whole-index change (follow-up 4), built from `git archive` in WSL
`~/grid-verify`. Database `db_grid`: that study's `make_variants.sh <protal> grid:20:` (its `a16` variant: 16 alleles,
the allele half's 20 genomes per focal species, no indexed alleles; now with grid levels). Runs: its `run_samples.sh`
with `g0:grid:--sparse_first,0`, `s4:grid:--sparse_first,4,--sparse_first_pick,stride`,
`g4:grid:--sparse_first,4,--sparse_first_pick,grid` (every sample grid-picked: the sparse index) and
`gw4:grid:--sparse_first,4,--sparse_first_pick,grid,--sparse_index,off` (the grid gate on the whole index); scored by
its `evaluate.py` into `~/strainfn/eval_grid*`.

| run | pairs given up | candidate alignments | records | calls: all / hard | 1-2 pairs | FP shallow / deep |
|---|---|---|---|---|---|---|
| no gate (`g0`) | 0 | 452,927 | 88,597 | 0.900 / 0.842 | 0.783 | 2 / 12 |
| stride C 4 (`s4`) | 124 of 88,800 | 452,764 | 88,592 | 0.900 / 0.842 | 0.783 | 2 / 12 |
| grid C 4, whole index (`gw4`) | 131 | 452,749 | 88,597 | 0.900 / 0.842 | 0.783 | 2 / 12 |
| grid C 4, sparse index (`g4`) | 131 | 557,934 | 88,597 | 0.867 / 0.743 | 0.675 | 1 / 7 |

- **No gate misses a strain.** The reads come from marker genes only, so the gates give up 0.15% of the pairs; the
  stride and the grid gate on the whole index call exactly what the ungated run calls.
- **The sparse index keeps the reads.** At the read level (deep samples) `g4` puts the hard strains' pairs on the
  right species more often than ungated (0.907 against 0.897; counted 0.807 against 0.798) and fewer on another
  species (0.032 against 0.038); easy and medium are the same within 0.003.
- **But the shared model calls fewer strains on it**, 0.867 against 0.900 (hard 0.743 against 0.842), with fewer
  false positives (1 against 2; deep 7 against 12). The 20 extra misses have counted pairs the model turned down ("1
  pair counted" 30 against 20, "2+ pairs counted, model said no" 27 against 17). At a lower knob `g4` is ahead at the
  same false positives (knob 0.2: 0.912, hard 0.851, 2 FP; `g0` at 0.5: 0.900, 0.842, 2 FP), so the evidence is there
  and the model's probabilities are lower. The single-seed anchors add 23% candidate alignments (456,986 single-seed
  anchors) for the same records, so the failed-candidate features rise: a model trained on ungated runs reads them as
  evidence against. A grid-picked database needs models trained on grid-picked runs, as `build_gtdb_database.py
  --sparse-first-pick grid` trains them (the collections seed as the sparse index, follow-up 4), and the r226 test
  rows (`--test-grid`) are scored by models trained otherwise, which will understate the grid's F1.
