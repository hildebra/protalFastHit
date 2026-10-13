# A sparse index by the genes' grid: a proposal

Date: 2026-10-13. Branch `audit-fixes`, on `0ab85f4` plus the uncommitted stride default
([r226 v27 gates](2026-10-13-r226-v27-c-gates.md), follow-up section). No code and no run: a design, with estimates
from earlier reports.

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
