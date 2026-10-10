# Ancestry sites and column weights on real GTDB r226 marker genes

**Data.**
- GTDB r226's marker-gene archives, GTDB-Tk's per-marker alignments of the representatives' proteins and GTDB's
  trees, downloaded unmodified to WSL `~/GTDB/r226/` (the server's layout; `MD5SUM.txt` with the server's sums).
  They stay packed; the chosen genomes' records are streamed out.
- [select_taxa.py](select_taxa.py): families with 4 or more genera of 4-12 species, two of them with 3 or more
  genomes (CheckM2 >= 90% complete, <= 5% contaminated); 5 families and 4 genera each (seed 1), every species of
  those genera, up to 20 strains per species.
  - Archaea (`ar53`): one family qualifies, Micrarchaeaceae (DPANN; 4 genera, 23 species, 54 strains, 53 markers).
  - Bacteria (`bac120`): 51 qualify; chosen Cyclobacteriaceae, Chitinophagaceae, Ruminococcaceae, Azotimanducaceae,
    Acidimicrobiaceae (20 genera, 147 species, 306 strains, 120 markers).
- [extract_subset.py](extract_subset.py): the chosen genomes' marker genes, their alignment rows and the tree, laid
  out as the server's extracted archives, in `~/real_ancestry/<set>`; [build_db.sh](build_db.sh) converts them
  (`gtdb_to_protal_db.py`) and builds a protal database with `--column_weights_alignment` `msa`, `aa` or `nt`;
  [validate.sh](validate.sh) / [validate_alignments.py](validate_alignments.py) compare each table's mapping with
  GTDB's alignment; [reference_choice.py](reference_choice.py) compares references for the alignments;
  [real_ancestry.py](real_ancestry.py) measures the ancestry sites through `scripts/ancestry_oracle.py` (the port of
  `AncestrySites.h`); [run.sh](run.sh) runs the selection, extraction and measurement;
  [placement_check.py](placement_check.py) checked the residue placement.
- protal: `0d961e1` (pooled codes), `c60c272` (protein alignments), then the GTDB-columns mode (this report's
  section 4), built in WSL `~/bidx` on 4 cores. The bacterial `*_validation.pairs.tsv` (9 MB each) are not kept;
  `validate.sh bac120 msa aa nt` writes them again.

**Questions (the user's, 2026-10-10).** Do the column weights' alignments work on real genes, and could they adapt to
divergent families? Does GTDB have a phylogeny to reconstruct genus and family ancestors from? How do the ancestry
sites fall on real genomes?

## Summary

- **The bases' alignments lose most of a real family; the proteins' half of it; GTDB's own columns nothing.** Within
  a family, the share of GTDB's aligned residue pairs that the table maps (coverage) and puts on the same column
  (agreement):

  | | archaea: coverage | agreement | bacteria: coverage | agreement |
  |---|---|---|---|---|
  | bases (`nt`, until `c60c272`) | 0.083 | 0.987 | 0.186 | 0.988 |
  | proteins (`aa`, `c60c272`) | 0.559 | 0.988 | 0.687 | 0.983 |
  | GTDB's columns (`msa`, now the default where the database has them) | 0.999 | 1.000 | 0.997 | 1.000 |

- **A reference central to its genus or family keeps every alignment short.** The share of species beyond 0.3
  protein divergence from their genus reference: random member 25.5% (archaea) / 11.6% (bacteria), medoid 12.4% /
  8.1%, the genus's ancestral sequence (Fitch on GTDB's tree) 5.9% / 2.7%; of genus references beyond 0.5 from the
  family's: 18.8% / 12.2%, 5.6% / 9.0%, 0% / 0.1%.
- **GTDB has the phylogeny:** `bac120_r226.tree` and `ar53_r226.tree`, maximum-likelihood trees of all
  representatives with the GTDB taxa on their nodes, and GTDB-Tk's alignments of every marker's proteins. The build
  now takes the columns from the alignment and each genus's and family's ancestral sequence from the tree (section
  4).
- **On real genomes the ancestry sites separate strains from congeners (AUC 0.994), the fixed sites perfectly
  (1.000), but a quarter of real strains carry the species' base at under 88.5% of its derived sites**: the sites
  come from one representative genome, so its own strain-level mutations count as the species'. The fixed sites
  (shared by the species' other known strains) remove those (section 5).

## 1. The column weights' alignments on real families

Archaea ([ar53_nt_build_lines.txt](ar53_nt_build_lines.txt), [ar53_aa_build_lines.txt](ar53_aa_build_lines.txt)):

| | bases (`nt`) | amino acids (`aa`) |
|---|---|---|
| species against their genus reference: kept | 126 of 372 | 450 of 553 |
| failed: too divergent / the aligner gave up | 107 / 139 | 75 / 28 |
| divergence of those compared, quartiles | 0.15, 0.20, 0.25 | 0.10, 0.20, 0.30 |
| genus references against the family's: kept | 53 of 126 | 93 of 126 |
| species' copies never aligned (their genus reference failed) | 326 | 145 |
| gene copies mapped | 225 (20 species) | 589 (23 species) |

Bacteria ([bac120_nt_build_lines.txt](bac120_nt_build_lines.txt), [bac120_aa_build_lines.txt](bac120_aa_build_lines.txt)):
gene copies mapped 6,327 of 16,080 by the bases (5,593 of 11,323 alignments failed), 12,668 by the proteins (1,372
of 13,443 failed; the genus references against the family's 0.25 / 0.35 / 0.50 of their amino acids apart).

The bases' limits (0.2 within a genus, 0.4 among genera) are below the real divergence: the species' nearest
congener is 0.12 / 0.18 / 0.25 of the bases away in the archaeal family, 0.07 / 0.13 / 0.21 in the bacterial ones
(quartiles, `congener_gaps`), and the genus reference is a random member, often farther.

## 2. Against GTDB's alignment

[ar53_*_validation.summary.txt](ar53_msa_validation.summary.txt), [bac120_*_validation.summary.txt](bac120_msa_validation.summary.txt):
every pair of representative copies of one family, GTDB's aligned residue pairs, by the pair's protein divergence on
GTDB's columns (archaea; coverage / agreement):

| divergence | pairs | `nt` | `aa` | `msa` |
|---|---|---|---|---|
| < 0.10 | 438 | 0.40 / 0.999 | 0.81 / 0.998 | 0.999 / 1.000 |
| 0.10-0.30 | 2,107 | 0.12 / 0.992 | 0.83 / 0.995 | 0.999 / 1.000 |
| 0.30-0.40 | 2,564 | 0.08 / 0.978 | 0.63 / 0.985 | 0.999 / 1.000 |
| 0.40-0.50 | 1,986 | 0.04 / 0.969 | 0.45 / 0.978 | 0.999 / 1.000 |
| 0.50-0.60 | 1,153 | 0.00 / 0.947 | 0.06 / 0.961 | 0.999 / 1.000 |

The protein alignments are as accurate as the bases' where both align and reach far more (accuracy slips above 0.4
of the amino acids: one mismatch penalty for all pairs, no substitution matrix). The `msa` agreement holds by
construction; it confirms that each copy's residues land on its own row's columns.

## 3. Which reference

[ar53_reference_choice.txt](ar53_reference_choice.txt), [bac120_reference_choice.txt](bac120_reference_choice.txt):
the protein divergence (GTDB's columns) of each species to its genus's reference, and of the genus references
(medoids) to the family's:

| | archaea: species, median / beyond 0.3 | genera, median / beyond 0.5 | bacteria: species | genera |
|---|---|---|---|---|
| random member | 0.212 / 25.5% | 0.406 / 18.8% | 0.119 / 11.6% | 0.300 / 12.2% |
| medoid | 0.135 / 12.4% | 0.317 / 5.6% | 0.088 / 8.1% | 0.292 / 9.0% |
| ancestral sequence (Fitch) | 0.111 / 5.9% | 0.164 / 0% | 0.075 / 2.7% | 0.111 / 0.1% |

The ancestor is reconstructed from the members it is compared with, so it is central by construction (that is
why it makes a good reference; the figures flatter it a little).

## 4. The column weights from GTDB's columns (`--column_weights_alignment msa`)

`ColumnWeightsMsa.h`; `gtdb_to_protal_db.py` writes `column_msa/<geneid>.faa.zst` (each representative's row of
GTDB-Tk's alignment, under its taxid) and `species_tree.nwk` (both trees, leaves renamed to taxids) from the release
when it has them (`download_gtdb.py` fetches them unless `--no_msa`); `--build` takes them by default (`auto`):
- each family's columns are the marker's alignment columns, three bases per amino-acid column; each copy's residues
  are placed on its row in order (`PlaceResidues`: a 4-residue window within 64 residues, an 8-residue one beyond, a
  residue just before an insertion continuing the previous one; a 4-residue window far ahead was a chance hit often
  enough to cost 1.4% of the archaeal pairs);
- within: every species of each genus votes (no genus reference); among: every genus of the family through its
  ancestral sequence (Fitch on the tree pruned to the genus); aa alike; the consensus base is the family's ancestral
  base where three genera or more carry a base.

[ar53_msa_build_lines.txt](ar53_msa_build_lines.txt), [bac120_msa_build_lines.txt](bac120_msa_build_lines.txt):

| | archaea | bacteria |
|---|---|---|
| copies mapped | 870 of 870 | 16,080 of 16,080 |
| the residues on a column (the rest GTDB's trimmed insertions) | 0.86 | 0.86 |
| genus ancestors, by parsimony on the tree | 167 of 172 | 2,364 of 2,378 |
| the family's ancestral base ambiguous at the root | 22% of columns | 22% |
| columns polarised | 31,974 | 598,440 |
| largest code within / among / aa | 6 / 4 / 7 | 7 / 4 / 7 |

No column reaches "conserved" (code 9) in these subsets: with the pooled pseudocounts that needs about 70 agreeing
species, and the subsets hold 4 genera of each family. At r226 every genus and species of a family votes.

## 5. The ancestry sites on real genomes

[ar53_ancestry.summary.txt](ar53_ancestry.summary.txt), [bac120_ancestry.summary.txt](bac120_ancestry.summary.txt). Per
species and marker, the species' sites from its representative against its congeners' representatives (protal's
consensus rule, without column weights); against them, each real strain of the species and each congener left out
of the sites in turn (a novel species the database lacks). Pairs with 10 or more sites covered:

| | bacteria: strains (306) | congeners (1,108) | archaea: strains (54) | congeners (146) |
|---|---|---|---|---|
| identity to the representative, median | 0.998 | 0.792 | 0.998 | 0.746 |
| agreement (the species' base), median (5%-95%) | 0.989 (0.678-1.000) | 0.204 (0.033-0.489) | 0.993 (0.984-1.000) | 0.277 (0.129-0.597) |
| the congeners' base, median (95%) | 0.008 (0.312) | 0.667 (0.942) | 0.004 (0.011) | 0.520 (0.765) |
| agreement at the fixed sites, median (5%) | 0.999 (0.909) | 0.023 (0.000) | 1.000 (0.997) | 0.065 (0.000) |
| AUC, strains over congeners: agreement / fixed sites | 0.994 / 1.000 | | 1.000 / 1.000 | |

By identity to the representative (bacteria; agreement of strains | congeners): 0.95-0.97 0.858 (16) | 0.218 (40);
0.97-0.99 0.851 (41) | 0.358 (4); >= 0.99 0.992 (249) | 0.989 (2).

- **Real strains 1-5% from their representative carry the species' base at only ~85% of its sites.** The sites are
  where the representative differs from its congeners, which includes the representative's own recent mutations:
  differences of one genome, not of the species. A strain that diverged before them carries the congeners' base
  there. The fixed sites (no other known strain of the species differs) are 0.999 for strains and 0.02 for
  congeners. This is why the plain agreement helped real strains nothing in r226 v22 (in-silico strains, mutated
  copies of the representative, keep all its sites) and why the fixed-site features carry the gain.
- **Sites per read** (bacteria, per species): 33.7 per kb compared (median; 5.4 at the 5th percentile), 5.3 per
  150-base read; 12% of reads cover no site (median species; a quarter of species more than 24%).
- **The hard case is rare here:** two congener pairs of the bacterial subset are >= 99% identical (agreement 0.989,
  indistinguishable); none of the archaeal ones is closer than 90%.

## 6. What to do next

1. **Build r226 with GTDB's columns** (`download_gtdb.py` fetches the alignments and trees; the converter writes them;
   `--build` takes them): the "Column weights codes" line shows whether conserved columns exist with whole families.
2. **The ancestry sites from more than one genome:** define the species' sites from the species' alleles (the strains
   in the full reference) as well as its representative, or weigh each site by how many of the species' genomes share
   it; the representative's private mutations then stop counting against the species' other strains.
3. **The ancestral genus sequence as the ancestry sites' reference** for species of small genera, where the
   congener consensus has too few votes.
