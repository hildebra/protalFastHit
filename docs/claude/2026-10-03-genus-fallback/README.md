# Reporting relatives of missing species at genus level: does it pay?

2026-10-03. Item 2 of the [v5/v6 report](../2026-10-03-r226-v5-v6-training/README.md)'s error budget: 74-88% of the
false positives at the knob curve are taxa of a genus whose species in the sample the database lacks. If protal reported
such a call as "an unknown species of genus G" instead of the congener, it would lose a false positive and keep a right
genus-level call; the ceiling was +0.009 (pe) species F1 for the 88 congener false positives of the paired-end test set.

Data: the r226 v5 run (`27423c6`, `+distance`, `local/v5`, git-ignored): its training and test tables, the
species-held-out scores of the training rows and the final model's test scores. `genus_fallback.py` makes
`genus_fallback_output.txt` (`python3 genus_fallback.py local/v5 "" _se _pb _ont`, from the repository root; ~5 minutes).

## Method

- Calls: each model at its knob curve, as protal calls by default.
- On the called taxa only, a second forest (300 trees, balanced classes, leaves of 3 or more) scores whether a taxon is
  *absent and of a genus with a species in the sample that the database lacks* (`meta_novel_congener`), from the
  normalised, adjacency and all relatives features, the first model's score and the sample's log10 fragments (all
  available to protal at run time).
- Converted: a called taxon at or above a threshold. It leaves the species calls; each sample gets one "unknown species
  of G" call per converted genus, right if the sample holds a species of G the database lacks.
- The threshold is chosen on the training rows by out-of-fold scores (5 folds grouped by species), maximising species
  F1; the forest is then refitted on all called training rows and applied once to the test set.

## Result: it does not pay

| read type | called test taxa | of them absent, genus with a missing species | threshold (training) | converted on the test set | species F1 |
|---|---|---|---|---|---|
| pe | 4,583 | 88 | 0.70 | 1 (genus right) | 0.9642 → 0.9644 |
| se | 4,399 | 67 | 0.75 | 0 | 0.9593 → 0.9593 |
| pb | 1,707 | 22 | 0.80 | 0 | 0.9761 → 0.9761 |
| ont | 1,735 | 25 | 0.70 | 0 | 0.9696 → 0.9696 |

On the training rows (out of fold), too, the best threshold converts nothing (pe, se, pb) or one taxon (ont):
**0-1 species per run would become genus calls, and species F1 changes by at most +0.0002.**

Lower thresholds, scored on the test set (not choices, the trade-off), pe:

| threshold | converted | absent, genus right | absent, genus wrong | present (species lost) | species F1 |
|---|---|---|---|---|---|
| 0.55 | 6 | 4 | 1 | 1 | 0.9647 |
| 0.45 | 16 | 8 | 1 | 7 | 0.9644 |
| 0.30 | 63 | 23 | 2 | 38 | 0.9626 |
| 0.20 | 146 | 40 | 5 | 101 | 0.9575 |

Every threshold that converts more than a handful converts more present species than false ones; se, pb and ont
show the same (`genus_fallback_output.txt`). The best test-set threshold (0.55, pe) gains 0.0005, chosen on the test set.

## Why: these false positives look like thin strains

The second forest ranks them well in the usual sense (AUC 0.92-0.95 out of fold) but they are 1.3-1.7% of the called taxa,
and AP is 0.19-0.26: no threshold catches many of them at a precision above one half. Medians of the called training
rows (pe):

| group | rows | fragments | identity | top identity | excess_median | genus_spill | em_own_share | first model's p |
|---|---|---|---|---|---|---|---|---|
| absent, genus has a missing species | 176 | 3 | 0.974 | 0.996 | 0.012 | 0.84 | 1.00 | 0.79 |
| absent, other | 196 | 1 | 0.990 | 1.000 | 0.002 | 0.48 | 1.00 | 0.90 |
| present strain, 1-10 fragments | 1,975 | 3 | 0.987 | 1.000 | 0.001 | 0.85 | 1.00 | 0.99 |
| present, representative | 5,313 | 34 | 0.990 | 1.000 | -0.001 | 1.80 | 1.00 | 1.00 |

- The reads of a missing species that reach a congener and pass the filters are its closest ones, mostly on
  conserved genes: 2-4 fragments at 0.97-0.98 identity, the best read 0.99+. That is what a present strain with 3
  fragments looks like (0.987, 1.000), and those strains are 11 times as many.
- `em_own_share` is 1 in every group: the missing species' reads have no other candidate in the database, so the read EM
  has nothing to split. `genus_spill` is the same as a thin strain's.
- The one signal is `excess_median` (reads differ more than their qualities explain; the second forest's top feature
  after `p`), which the first model already uses. Where a missing species leaves more reads (10 or more fragments) the
  first model rejects it already: 89-90% of the remaining false positives have fewer than 10.

So a genus-level fallback chosen from these features only re-draws the first model's threshold among thin taxa. The
error-budget ceiling (+0.009 pe) assumed a perfect separation the features do not give.

## What could still separate them

- **More reads per decision:** the missing species' reads are spread over its genes but only the conserved ones reach
  the congener. A feature of *where* on the gene the reads lie (the variable regions of a congener gene have none) or of
  the identity of the reads on fast genes alone might separate them; `relative_close_share` tries this by gene pairs and
  helps little here.
- **The strain genomes of item 1**: with strains in the database, a present strain's reads would be ~0.99+ identical
  and no longer look like a missing congener's, which would also make these false positives separable.
- Not a rule on the current features; no protal change suggested.

The website is not affected.
