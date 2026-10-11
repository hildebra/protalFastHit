# NCBI taxids for GTDB species: what GTDB provides, a majority table in the database, and placeholders

2026-10-11. Code base `f4850d1` (branch `audit-fixes`) plus the change described here, uncommitted
when this was written. Data: GTDB r226 metadata (`bac120_metadata_r226.tsv.gz`,
`ar53_metadata_r226.tsv.gz` in WSL `~/GTDB/r226`): 732,475 genomes, 143,614 species. Machine:
the WSL box, 4 cores (`taskset -c 0-3`).

## 1. What GTDB provides

GTDB has no table from its species to NCBI taxids, and its species have no numeric ids, only
names. Its metadata gives every genome NCBI's view of it:

| Column | |
|---|---|
| `ncbi_taxid` | the genome's own taxid at NCBI (strain level where NCBI has one) |
| `ncbi_species_taxid` | the taxid of the NCBI species it belongs to; filled for all 732,475 r226 genomes |
| `ncbi_taxonomy_unfiltered` | NCBI's lineage as it is, with its placeholder names ("s__uncultured Odoribacter sp.") |
| `ncbi_taxonomy` | GTDB's filtered copy of it, with every name that is no proper taxon left empty ("s__") |
| `ncbi_organism_name` | the organism name of the assembly |

The species translation has to be made by grouping genomes by GTDB species. GTDB-Tk's
`gtdb_to_ncbi_majority_vote.py` does that for classified genomes. `gtdb_to_taxdump` (Youngblut)
makes NCBI-style `names.dmp`/`nodes.dmp` with made-up ids for GTDB taxa.

## 2. The majority table

The converter (`scripts/mini_db/gtdb_to_protal_db.py`) now writes `ncbi_taxids.tsv`. For each
species it takes the NCBI species taxid that most of the species' genomes have. A tie goes to the
representative's taxid, else to the smallest. `named` is 1 where GTDB's filtered NCBI lineage keeps
the species name, 0 for an NCBI placeholder. The other columns are `ncbi_name` (from the unfiltered
lineage), `ncbi_genomes`/`genomes` (the votes and all the species' genomes) and `shared_by` (the
species of the release whose majority is the same taxid).

r226, by the converter's own functions (`read_ncbi_species`, `ncbi_majority`) on the metadata,
with species and taxids as a conversion gives them (`r226_table.py`, section 7):

| Majority taxid | Species | |
|---|---|---|
| named NCBI species, this species' alone | 17,617 | 12.3% |
| named NCBI species, shared with other GTDB species | 7,547 | 5.3% |
| NCBI placeholder, this species' alone ("Vibrio sp. HB236076") | 16,440 | 11.4% |
| NCBI placeholder, shared ("bacterium", "uncultured Lachnospiraceae bacterium") | 102,010 | 71.0% |

- 82.5% of the species get a placeholder; they hold 44% of the genomes (410,745 of 732,475 genomes
  are in species with a named NCBI species).
- How many species share a taxid: 1 for 34,057 species (23.7%), 2 for 4,172, 3-10 for 11,100,
  11-100 for 29,132, 101-1,000 for 42,142, more than 1,000 for 23,011 (16.0%). The most shared are
  "uncultured Clostridia bacterium" (244328, 2,311 species), "bacterium" (1869227, 2,087) and
  "Pseudomonadota bacterium" (1977087, 2,078).
- Named species are lumped where GTDB split them: *Escherichia coli* (562) is the majority of
  eight GTDB species (*E. coli*, *E. coli_E*, *E. coli_F* and five MAG species), *Bacillus cereus*
  (1396) of 19, *Akkermansia muciniphila* (239935) of 58.
- Support: 94,012 species have one genome (their label is the genome's). Of the 49,602 with more,
  19,137 agree on one taxid, 30,263 have a majority above half, 11,917 a half and 7,422 below
  half. 14,895 species had a tie, broken by the representative's taxid.
- The table is 12.8 MB at r226 (1.9 MB at zstd level 19), read only by runs with `--ncbi_taxids`.
  Reading the metadata and voting add ~9 s to a conversion (6.5 s and 2.3 s).

### Mislabels at NCBI

The vote takes NCBI's labels as they are. Bulk MAG submissions are sometimes labelled with another
organism. Of the 58 GTDB species with "Akkermansia muciniphila" as their majority, most are MAGs
of other genera and families: *Prevotella* sp963621665, *Succinivibrio* sp963622535, *Rhodoferax_C*
sp943327015, *Luteolibacter muciniphilus_A* and others; nine are *Akkermansia* species. Most of the
others rest on one genome, or on a tie that the representative's label won. The archaeon
*Methanarcanum* sp963632775 gets "Bacteroides salyersiae" (1 of 9 genomes, a tie).

Of the 25,164 named picks, 1,326 (5.3%) have both an NCBI genus and an epithet other than the GTDB
species' (941 of them from one genome). That counts spelling and gender variants too
(*cuticulare*/*cuticularis*), so it is an upper bound. NCBI's and GTDB's lineages differ too often
in their class and order names to detect mislabels by comparing lineages: they part at the domain
for only 5 named picks.

### A "named first" vote, tried and dropped

Before choosing the plain majority I tried a vote among the genomes with a named NCBI species,
falling back to all genomes (scratch `ncbi_analysis.py`; ties there to the representative, then
the smallest id). It changed 3,539 species, where a named minority beat a placeholder majority.
Its NCBI genus differed from the GTDB genus in 1,088 of them (31%), e.g. *PeH17* sp000435055 →
*Fusicatenibacter saccharivorans* (1 of 112 genomes), *Limnohabitans* sp945870575 → *Bordetella
parapertussis* (1 of 6). A minority named label is as often a mislabel as a better name. The
plain majority stays, as asked.

## 3. Placeholders for species NCBI does not name

Every GTDB species has a taxid in the table: no genome lacks an NCBI species taxid. "Not covered"
means that the taxid is one of NCBI's placeholders. The question is what the profiles should
print for them:

| Option | NCBI id? | Unique to the GTDB species? | Information |
|---|---|---|---|
| a. NCBI's placeholder taxid (implemented) | yes, NCBI's own for these genomes | for 16,440 of 118,450 placeholder species | NCBI's tools place it under the deepest rank NCBI names (below) |
| b. the deepest named NCBI ancestor's taxid (as MetaPhlAn 4 does for its unknown SGBs) | yes | less often than a | the same rank as a, without a fake species; needs NCBI's taxdump at build time to turn names into ids (with homonyms resolved), at a newer date than GTDB's snapshot |
| c. made-up ids (negative, or an offset range) | no | yes | none for NCBI's tools |
| d. the GTDB name | no | yes | full GTDB lineage; already the profile's second column |

The deepest rank NCBI names for the placeholder picks (GTDB's filtered NCBI lineage of the
genomes with the taxid): genus 32,339 (27.3%), family 23,745 (20.0%), order 16,361 (13.8%), class
17,600 (14.9%), phylum 23,026 (19.4%), domain 5,379 (4.5%).

Option a is better than the GTDB name where an NCBI id is wanted. It is a valid id, it is what
NCBI itself says about these genomes, and its lineage resolves in NCBI's tools (`taxonkit`,
Krona). For 16,440 species it is even the species' own. The cost is lumping: a table keyed by
NCBI taxid merges up to 2,311 GTDB species into one row. The profile line keeps the GTDB lineage
in its second column, so no row is lost. The NCBI name column makes a placeholder visible
("uncultured Clostridia bacterium"), and `named` and `shared_by` in the database's table say it
for scripts. Option b gives the same information with fewer, coarser ids, at the price of a taxdump
download and name resolution in the build; it is worth doing only if downstream tools must not see
placeholder species.

## 4. In protal

- `ncbi_taxids.tsv` (`src/Taxonomy/NcbiTaxids.h`): read by its header (`taxid`, `ncbi_taxid`,
  `ncbi_name`, `named`; the other columns skipped). `--build` packs it, `--unpack_db` writes it,
  and `--add_tables` stores it after checking every taxid and representative against the
  database's taxonomy, as for `species_priors.tsv`. `gtdb_to_protal_db.py --priors_only` writes it
  beside the priors, so a database converted before today gets it without a rebuild. A
  `--from_db` copy (the GTDB build's training database) copies it.
- `--ncbi_taxids`: each `.profile` line ends with the taxid and NCBI's name; `NA NA` for a species
  the table lacks; `? ?` on the unknown-share line. The other columns are unchanged. A database
  without the table stops the run with exit 2 before anything is aligned; a broken table exits 8.
  Without the option the table is not read and the profiles are byte-identical.
- The synthetic GTDB release (`simulate_gtdb_release.py`) gets made-up NCBI columns:
  every third species, and the other species' genomes from the third on, have an "uncultured
  <genus> sp." placeholder of their genus. The other random draws are unchanged.

Tests: C++ unit tests 501 passed (new: `NcbiTaxids.*`, 3), `test_mini_db.py` 21 passed (new:
`test_ncbi_taxids`, `test_ncbi_majority`), `test_gtdb_build.py` 27 passed, e2e 150 passed
(new: `ReadTypeModelTest.test_ncbi_taxids`: the columns on the mini database, a database without
the table, `--add_tables` of another release's table refused and of a renamed one used),
`test_gtdb_pipeline.py` 7 passed (section 6).

## 5. Not done

- No check of the picks against GTDB's lineage (section 2: too noisy above the domain); `ncbi_genomes`
  and `genomes` are there to filter on.
- Option b of section 3.
- NCBI taxids for higher ranks (a CAMI-format profile with a TAXPATH) would need the taxdump too.
- The website does not list `--ncbi_taxids` or the profile's extra columns.

## 6. Pipeline test

`test_gtdb_pipeline.py` (the GTDB build end to end on the synthetic release, whose metadata now has
the NCBI columns, so its converter writes the table and both builds pack it): 7 tests passed in
228 s.

## 7. Commands

```bash
# build and tests (WSL; a copy of the checkout in ~/protal-ncbi/src)
cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DPROTAL_BUILD_TESTS=ON
taskset -c 0-3 nice -n 5 ninja -C build -j4 && (cd build && taskset -c 0-3 ctest -j4)
taskset -c 0-3 python3 -m unittest scripts/mini_db/test_mini_db.py
PROTAL=build/protal scripts/mini_db/build_mini_db.sh mini_db
PROTAL_TEST_DB=mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
    python3 -m unittest tests/e2e/test_protal_e2e.py
PROTAL=build/protal SIMULATE=build/simulate_metagenomes PROTAL_TRAIN_PYTHON=~/micromamba/envs/protal-db-build/bin/python \
    PROTAL_TESTS_REQUIRED=1 python3 -m unittest scripts/mini_db/test_gtdb_build.py scripts/mini_db/test_gtdb_pipeline.py
```

The r226 numbers came from scratch scripts that import the converter's functions (the
representatives' marker files, which a real `--priors_only` run reads to know the species, are
packed on this machine): `r226_table.py` (the table and section 2's counts), `r226_suspect.py`
(the mislabel bound), `r226_consistency.py` (lineage agreement), `r226_depth.py` (section 3's
ranks), `ncbi_analysis.py` and `ncbi_named_first.py` (the first look and the "named first" vote).
They were not kept.
