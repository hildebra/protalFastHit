# Congener gaps, foreign rates and the reads lost before alignment: implementation

**Question (2026-10-07).** A follow-up to [the reads behind the errors](../2026-10-07-error-read-signatures/README.md):

- Align each marker gene within its genus at build time, keep each copy's min and median distance to its congeners, and
  use them as thresholds for counting a taxon's hits, as features.
- Pass a per-copy foreign-read rate, learned from a scan of the genomes, as another feature.
- Recover the missed strains' reads that never reached the aligner, which the analysis traced to `--align_top` 3 in
  crowded genera.

The user chose:
- implement the gaps directly;
- take the foreign rates from a tiled self-alignment of the genomes at hand;
- add adaptive candidates and an untried-candidate feature, rather than an order-dependent prior;
- build and test locally.

**Short answer.** All of it is implemented, tested here and in the default feature set, but not yet trained at GTDB
scale: the next r226 build is the first test of its worth. Branch `congener-gaps` (from `0aea3fb`).

## What it does

### 1. The gaps to the congeners' copies (`congener_gaps.tsv`, `gaps` features)

**At build time.** `--build` (`WriteCongenerGaps`, `src/SequenceUtils/CongenerGaps.h`) goes gene by gene. Each
species' copy is aligned against the copies of its genus's other species:
- all of them up to 24 others;
- in larger genera, its 4 nearest by the 128-hash k-mer sketches (for the minimum) and 16 others drawn by a hash of the
  pair (whose median gives the median).

The alignments use WFA2 with protal's read scores (4/6/2), ends free up to 15% of each copy (congeners' gene calls
start and end in different places). Both copies must overlap over half the shorter one, and an alignment beyond a
distance of 0.4 is given up.

The distance is the share of aligned columns that differ (mismatches and gap bases), as a read's divergence is
counted. Per copy, the table keeps:
- the nearest congener's distance (min);
- the median of the sample's distances (median);
- the number of congeners with the gene.

They are stored in 1/10,000 units, one line per species. A run loads the table into the `GenomeLoader` in parallel
with the other tables: 8 bytes per copy, about 120 MB at r226.

**Features** (per taxon, from its kept records, MAPQ 4 or more; −1 without the table or without a qualifying record):

| feature | what |
|---|---|
| `gap_informative_share` | share of the kept records on copies whose nearest congener is at least 0.005 away (nearer, a read cannot tell them apart) |
| `gap_within_min_share` | of those, the share whose divergence is below the copy's min |
| `gap_within_median_share` | the share below the copy's median |
| `gap_position` | median of divergence over the min, from a 30-bin histogram of 0.1 (integer counts, so the same on any thread count) |

These are the user's "min and median as thresholds for counting hits", as shares. A strain's reads should stay within
the gap to the nearest congener. A novel congener's reads should lie about as far from the reference as its congeners
do. Unlike the run's `relative_*` features, the table knows every congener in the database, not only those present in
the sample.

### 2. The foreign rates (`foreign_rates.tsv`, `foreign` features)

**The scan.** `simulate_metagenomes --tiles LENGTH:STRIDE` cuts error-free reads from every contig of every genome of
the table:
- reads named `<genome>:<contig>:<start>`;
- the species of `--tile_exclude` left out;
- genomes read in parallel, written in table order, so the same file for any thread count.

[`scripts/foreign_rates.py`](../../../scripts/foreign_rates.py) then:
1. runs that (150 bases every 500);
2. aligns the tiles once with protal (single-end, `--no_profile`);
3. counts each read's best record (MAPQ 4 or more) for its copy: own species, another species of the genus, or another
   genus, by the genome table's lineage against the database's taxonomy;
4. writes `taxid<TAB>gene:reads:foreign:foreign_genus,...`.

`protal --add_tables foreign_rates.tsv --db DB` stores it, the same way `--add_model` stores models (members replaced or
added, the models last, in place where possible). `--add_tables` also takes `congener_gaps.tsv`.

**In the build.** `build_gtdb_database.py` runs the scan right after the training database is built:
- the held-out species are left out, so the table knows nothing of the species the models learn to find missing;
- the table is stored in the training database and copied into the finished database's folder before its build packs it
  (the two share taxids);
- `--no-foreign-rates` skips it, `--foreign-stride` sets the spacing.

**Features** (−1 without the table or without a record on a scanned copy):
- `foreign_scanned_share`: the share of kept records on scanned copies;
- `foreign_copy_share` and `foreign_genus_copy_share`: those copies' mean shares of other species' and other genera's
  reads, foreign / (reads + 1).

### 3. The reads lost before alignment (`ZC`, `--adaptive_candidates`, `untried` feature)

The two alternatives to an order-dependent prior that the user chose:

- **Untried candidates.** A short read's first record gets `ZC:Z:<taxid>,...`: up to 8 taxa of its crowd (anchors at
  least 0.8 as long as its longest, as for `ZN`) that it was never aligned against, the strongest first. A pair joins
  both mates' and leaves out what either tried or aligned to. They travel as flagged entries of the failed-candidate
  list: `ZF` and the header's counts of unaligned reads leave them out, and an unmapped record is written only for real
  failed candidates.
  - The profiler counts them per taxon that has records: `untried_candidate_rate` = untried / (untried +
    fragments_all), as `failed_candidate_rate` does.
- **Adaptive candidates.** A short read whose best alignment so far is divergent (score identity below 0.99) also tries
  more anchors beyond `--align_top`, up to `--adaptive_candidates` (default 7). They are taken from its crowd, in anchor
  order, of taxa in that best alignment's genus not tried yet.
  - This is per read, so the result does not depend on read order or threads.
  - The run loads the taxonomy for the genus map, also when it only aligns.
  - The log's `adaptive candidates:` line counts the extra alignments. `--adaptive_candidates 0` aligns as before.
  - Long reads are not changed.

## Tests (WSL, 4 cores, niced)

- C++ unit tests, `tests/test_CongenerGaps.cpp`, 8 new tests:
  - the alignment distance (substitutions, free ends, a gap, unrelated copies);
  - the scan of a genus, all pairs and sampled, the same on 1 and 3 threads;
  - both tables' write, read and errors;
  - the untried tags;
  - the handler: a read of a strain ranked fourth behind three congeners aligns to a congener and leaves its species in
    `ZC` with `--adaptive_candidates 0`, and finds its species with 7; a read of the reference tries nothing more;
  - the profiler's 8 features with and without the tables, on 1 and 3 threads.
- All 446 unit tests ran: 443 passed, 2 were skipped, and one failed, `IlluminaSimulation.AFailedStreamIsCutOff`. That
  test is a named-pipe timing test this change does not touch; it passed 3 runs of 3 on its own.
- `scripts/test_foreign_rates.py` (2 tests): the genome and taxonomy readers, the counting of a scan's SAM (best record,
  MAPQ, unmapped, unknown genome, supplementary first), and the table format.
- `scripts/test_model_pmml.py`: the default feature set pinned with the three groups, and a table of a protal before them
  failing clearly.
- The CI suites on this build, with `PROTAL_TESTS_REQUIRED=1`, all passed:
  - the mini database, whose build wrote `congener_gaps.tsv`: 222 copies of its 2 congeners, nearest distances
    0.075-0.086, in 28 ms;
  - the end-to-end tests (142);
  - the mock community (`examples/mini_db/run.sh`);
  - the mini database and GTDB build tests (76). They include the whole build with the scan step: `foreign_rates.log`,
    the stage, the table in the training database, and the training tables' `gap_informative_share`,
    `foreign_scanned_share` and `untried_candidate_rate` known (not −1) where reads landed.
  - the script and trainer tests (76).

## Cost at GTDB scale (estimates, not measured)

- **Congener gaps.** r226 has 143,614 species in 29,405 genera, the largest with 1,739 species. With the caps that is at
  most ~1.5×10⁸ gene alignments: genera of 25 species or fewer in full, larger ones at 20 per copy, ~101 genes per
  species. The mini database aligned 111 pairs at ~8% divergence in 28 ms on 4 threads. At that rate the step takes
  roughly 10 CPU-hours, about 12 minutes of a 52-thread build. The build log's `Congener gaps took` line will say.
- **Foreign rates.** The scan reads every genome at hand once: ~100 Gbp at r226 in 150-base tiles every 500 bases is
  ~2×10⁸ reads, one protal run. At r226's measured single-end speed that is minutes to tens of minutes. It runs between
  the training database's build and the profiling, so it lengthens the build's critical path by that much.
- **Adaptive candidates.** Extra alignments only for divergent reads in crowded genera. The handler's k-mer screen
  refuses 90% of candidates at r226 for ~0.7 µs each, so most of the cost is screening; the run's log counts them.
- **Memory at run time.** ~120 MB for the gaps and up to ~120 MB for the foreign rates (only the scanned copies). Each
  table is held as text while it is parsed: ~350 MB and less.

## Overlap with the ancestry sites (the other branch)

A parallel session implemented "ancestry sites" (`AncestrySites.h`, uncommitted in the main checkout on 2026-10-07).

**How they work.**
- A run compares a reference's copy of a gene with its nearest congener's copy (from `species_neighbours.tsv`, the first
  of 3 that pairs). It pairs the two copies' shared unique 12-mers on the main diagonal, without an alignment.
- It keeps the positions where they differ and the congener's base there. This is computed lazily per (species, gene)
  a sample touches and cached.
- Per read: the sites it covers, and whether it has the species' base or the congener's (`ancestry_sites_per_record`,
  `ancestry_agreement`, `ancestry_congener_share`).

**Where the two overlap.** Only part 1 of this branch, and only half of it:

| question | gaps (this branch) | ancestry (the other) |
|---|---|---|
| is the read on the species' side of its nearest congener? | `gap_within_min_share`, `gap_position`: the read's divergence over the whole gene against the copy's distance to the nearest congener | `ancestry_agreement`, `ancestry_congener_share`: base by base at the sites where the two differ |
| how far is a typical congener? | `gap_within_median_share` (16-24 congeners per copy) | no |
| can the gene tell at all? | `gap_informative_share` (nearest congener ≥ 0.005 away) | `ancestry_sites_per_record` (sites per record), close in spirit |
| which congener | the nearest by that gene's own sketches, any congener of the genus | the nearest by all markers, within 0.15, from species_neighbours (≤ 16) |
| computed | at build time, WFA2 alignment with indels counted; stored per copy (~120 MB) | at run time, per touched copy; k-mer pairing, stretches across indels not compared |

**Where it overlaps,** the ancestry sites ask the same question better: position by position, they avoid the caveat
this branch lists, that a gene's divergence varies along it. `gap_within_min_share` and `gap_position` are probably
redundant with `ancestry_agreement` once both are trained. The median and the informativeness are not covered by the
ancestry sites. Neither are parts 2 and 3: the foreign rates, `--add_tables`, `ZC`, the adaptive candidates and
`untried_candidate_rate`.

**Duplicated work.** The two branches compute the nearest-congener comparison twice, at build time here and at run time
there, and pick the nearest congener differently. Options:
1. Train both and let the next r226 build's ablation (varimp, `--features` without `gaps`) decide. That is the cheapest
   now, and correlated features cost a gradient-boosted model little.
2. Keep only the median and the informative share in `gaps`. Drop `gap_within_min_share` and `gap_position`, and with
   them the need to find each copy's nearest congener: the build would align only the 16-24 sampled congeners.
3. Give `congener_gaps.tsv` each copy's nearest congener by that gene, 4 more bytes per copy. The ancestry cache would
   then compare against the gene's own nearest congener instead of the species-level one, including congeners beyond
   species_neighbours' 0.15. That makes the two consistent and drops the guess of 3 tries.

The recommendation is 1 now and 2 or 3 after the build. Both branches touch the same files (RecordEvidence,
TaxonFeatures, `model_features.py`'s default set): merging them means keeping both sides.

## What is not done

- No r226 build: none of the 8 features has been trained or shown to help.
- `--align_top` 5 or 10 has not been timed.
- Long reads get neither `ZC` nor adaptive candidates.
- A read's divergence is compared with the whole gene's gap, though the gap varies along the gene. A per-window table
  (about 20 bytes per copy) would make the test sharper.
- The foreign rates count reads of every genome at hand, not only reads from outside the marker genes. A transferred
  gene's reads and a congener's conserved-gene reads are counted alike.

## Commands

    cmake --build build --target protal simulate_metagenomes protal_tests && build/tests/protal_tests --gtest_filter='CongenerGaps*:ForeignRates*:UntriedCandidates*:AdaptiveCandidates*:CopyFeatures*'
    python3 scripts/foreign_rates.py --db DB --genome-table genomes.tsv --taxonomy internal_taxonomy.dmp --exclude heldout_species.txt --out foreign_rates.tsv -t 32
    protal --add_tables foreign_rates.tsv --db DB
