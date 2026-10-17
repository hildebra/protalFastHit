# Does protal fit its training scenarios? Leave one scenario out, and communities never trained on

The user asked whether protal overfits its training scenarios. They asked for a cross-validation in the database build
that tests a very different community composition.

## 1. What the build checked before

| Check | What it holds out | Same kind of community as the training? |
|---|---|---|
| Cross-validation (`--cv-scheme samples+species`) | a row's sample and its species | yes |
| Clades held out (`--evaluation full`) | whole genera to phyla | yes |
| Independent test set (`--test-*`) | another design: depths, read setups, sizes | yes, the design's kind |
| Scenario hold-out samples (`--scenario-test-samples`) | the scenario's samples of another seed | yes, the same scenario |

None of them scores a kind of community the models never saw. A model that has learnt the five scenarios' communities
(gut, moderate, soil, shallow soil, host) rather than what they share looks as good as one that generalises.

## 2. What is added (both checks, on by default)

The user chose both checks, on by default, and two out-of-distribution kinds: dense congeners and extreme unevenness.

**Out-of-distribution scenarios.** `scripts/scenarios.py` gains `OOD_PRESETS`.
`build_gtdb_database.py --ood-scenarios` (default `congeners_dense,uneven_extreme`; `none` for none) puts them into
the test set alone: `--scenario-test-samples` hold-out samples each, no hold-in samples.

| Scenario | Community | Reads |
|---|---|---|
| `congeners_dense` | 100-400 species, nine in ten in groups of 4-12 of one genus (`congeners 0.9:4-12`), second and third strains in half and a quarter, 10% lacking | Illumina PE 150 Q35 at 10M pairs, Ultima at 10M reads, PacBio and Nanopore at 3 Gb |
| `uneven_extreme` | 500-3,000 species at lognormal sigma 3.0 and 3.5 (the training's reach 2.5), 30% lacking | Illumina at 20M pairs, Ultima at 20M reads, PacBio and Nanopore at 6 Gb |

- **Never trained on:** `--scenarios all` leaves them out, and naming one in both options is an error.
- **They choose nothing:** `--features auto` reports their F1 per candidate set but chooses by cross-validation alone;
  the knob and the model are chosen without them.

**Leave one scenario out.** `machine_learning_cmdline.py --scenario-transfer` (on; `--no-scenario-transfer`):

- **The fits:** for each scenario with hold-in rows, a model of the final one's settings and row weights fitted
  without them, and one fitted on the design's rows alone.
- **The scoring:** every scenario's hold-out samples are scored by these models beside the final model, at the final
  model's knob (or knob curve).
- **The output:**
  - a report section, "Scenarios: carried over to communities not trained on", with a table (F1, precision,
    sensitivity, FP per sample, FN rate, best F1) and a line per scenario;
  - `metrics.json` `scenario_transfer`;
  - `build_metadata.tsv` rows `model_<type>_scenario_transfer`, `classifier_ood_scenarios` and
    `classifier_scenario_transfer`.
- **In the build:** each read type's model of the database's gate and `--model`'s kind runs it; the gate models and
  the other kind pass `--no-scenario-transfer`. The build's `--no-scenario-transfer` turns it off everywhere.

**How to read it:**

- **A scenario trained on:** the final model's hold-out F1 against the model fitted without that scenario. The gap is
  what its own samples bought, and so how much the models depend on having seen that kind.
- **A scenario never trained on:** its F1 against the in-distribution scenarios' hold-out F1 says how the models do
  on something new.
- **The design-alone model:** shows what the scenarios as a whole add, on every kind.

A large gap means the models fit the kinds of communities they saw; a small one means they carry over.

## 3. Cost

- **Simulation:** two more scenarios' hold-out samples, 4 each per read type at 10-20M pairs (3-6 Gb of long reads):
  about two fifths more scenario test samples than the five scenarios' 20.
- **Training:** one fit per scenario and one more per read type, each the final model's fit. At r226 v18 a GBM fit
  took 4-26 s; the v31 tables are larger, so expect minutes per read type.

## 4. Checks run

- **`test_model_pmml.py` `TrainerScenariosTest.test_scenario_transfer`:**
  - gut (trained) gets the final model's row, the model's without gut's hold-in rows, and the design's;
  - host (test table alone) gets the final model's and the design's rows, marked never trained on;
  - the final model's row equals the scenario table's hold-out row;
  - `--no-scenario-transfer` drops the section.
- **`test_collector.py`:** the OOD presets are well-formed, selected by name, and never by `all`.
- **`test_gtdb_build.py`:** 46 pass.
- **`test_gtdb_pipeline.py` `test_f_scenarios`** (a small `--scenario-file`, the two OOD scenarios small too):
  - they are in the test set's paired-end table, and in neither training table;
  - the pe report has the section and the lines `scenario <name> (never trained on): ...`;
  - `metrics.json` and `build_metadata.tsv` have them;
  - the feature choice reports them without choosing by them.

  It passes (112 s). In the full pipeline run before (8 tests, 534 s) the other 7 tests passed. `test_f` failed only on
  the error reads' list of scenarios: it now names the out-of-distribution test samples too, and the test expects that.

The binaries were those of commit 830a92b (WSL `~/contam-build`); the scripts were the working tree's.

## 5. Not done, and caveats

- **No r226 build has run it:** the gaps are unmeasured.
- **What the out-of-distribution scenarios cover:** they differ in community composition only. They are still built
  from GTDB genomes with the same read models, so a technology protal never trained on, or genomes far from any in the
  release, are not tested.
- **Leave one scenario out measures transfer between the five kinds:** a gap there is a lower bound on what an unknown
  kind of real study costs.
- **The knob is the final model's:** a model without a scenario might choose another knob, so `best F1` (the
  threshold best for those rows) is given beside the F1 at the knob.
