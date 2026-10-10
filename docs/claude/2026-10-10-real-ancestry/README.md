# Ancestry sites and column weights on real GTDB r226 marker genes

**Status:** in progress (2026-10-10). The archaea are done; the bacterial subset waits for GTDB's archives.

**Data.**
- GTDB r226's marker-gene archives, downloaded unmodified to WSL `~/GTDB/r226/` (the server's layout, `MD5SUM.txt`
  with the server's sums). They stay packed; the chosen genomes' records are streamed out.
- [select_taxa.py](select_taxa.py): families with 4 or more genera of 4-12 species, two of them with 3 or more
  genomes (CheckM2 >= 90% complete, <= 5% contaminated); 5 families and 4 genera each (seed 1), every species of
  those genera, up to 20 strains per species. Archaea (`ar53`): one family qualifies. Bacteria (`bac120`): 51
  qualify; chosen Cyclobacteriaceae, Chitinophagaceae, Ruminococcaceae, Azotimanducaceae, Acidimicrobiaceae (20
  genera, 147 species, 306 strains).
- [extract_subset.py](extract_subset.py): the records of the chosen genomes, laid out as the server's extracted
  archives, in `~/real_ancestry/<set>`; [build_db.sh](build_db.sh) converts them (`gtdb_to_protal_db.py`) and builds
  a protal database with the strain alleles and the column weights; [real_ancestry.py](real_ancestry.py) measures
  the ancestry sites (through `scripts/ancestry_oracle.py`, the port of `AncestrySites.h`); [run.sh](run.sh) runs
  the selection, extraction and measurement.
- protal at `0d961e1` plus the amino-acid alignment of the column weights (built in WSL `~/bidx`, 4 cores).

## 1. The column weights' alignments on a real family

The archaeal family Micrarchaeaceae (DPANN; 4 genera, 23 species, 54 strains, 53 marker genes;
[ar53_nt_build_lines.txt](ar53_nt_build_lines.txt), [ar53_aa_build_lines.txt](ar53_aa_build_lines.txt)):

| | bases (`nt`) | amino acids (`aa`) |
|---|---|---|
| species against their genus reference: kept | 126 of 372 | 450 of 553 |
| failed: too divergent / the aligner gave up | 107 / 139 | 75 / 28 |
| divergence of those compared, quartiles | 0.15, 0.20, 0.25 | 0.10, 0.20, 0.30 |
| genus references against the family's: kept | 53 of 126 | 93 of 126 |
| divergence of those compared, quartiles | 0.40, 0.45, 0.50 | 0.40, 0.45, 0.55 |
| species' copies never aligned (their genus reference failed) | 326 | 145 |
| gene copies mapped | 225 (20 species) | 589 (23 species) |
| within codes: columns with an estimate, largest | 26,397, 5 | 34,623, 6 |

- The bases' limits (0.2 within a genus, 0.4 among genera) are below this family's real divergence: the species'
  nearest congener is 0.12 / 0.18 / 0.25 of the bases away (quartiles, `congener_gaps`). Two thirds of the species'
  alignments failed, and every species of a genus whose reference failed was lost with it.
- As proteins most copies align. The family is extreme even so (its genus references differ from the family's by
  0.40-0.55 of their amino acids), near where one mismatch penalty for all amino-acid pairs aligns poorly.
- No column reaches "conserved" (code 9): with the pooled pseudocounts that needs about 70 agreeing species, and
  the family has 23.

## 2. The ancestry sites on real genomes (archaea)

[ar53_ancestry.summary.txt](ar53_ancestry.summary.txt): at the species' derived sites (protal's consensus rule,
congeners as compared by `CompareCopies`), the species' strains carry the species' base at 0.98-1.00 of the sites
(median 0.993), its congeners, each left out of the sites as a novel species would be, at 0.13-0.60 (median 0.28);
the AUC is 1.0. But this family's strains are >= 99% identical to their representatives and its congeners about
75%: the hard case of protal's errors, a novel congener 95-99% identical, does not occur in it.
