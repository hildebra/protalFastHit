#!/usr/bin/env python3
"""database_check.py - a test set's samples profiled again, on another database or with another protal option, scored
beside the calls the trainer made on the same samples (PREFIX.calls.tsv.gz, its test rows: the final model's calls as
protal makes them by default, on the training database). build_gtdb_database.py runs one such check, the final check
(model_logs/final_db_check.tsv): the test samples written to the disk profiled again on the finished database, which has
every species and alleles from every genome, with its default models at their gates. Only the rows of species that
were not held out count (neither held-out species nor those of held-out clades: the training database lacked them, and
so they have no row there), on both databases. (The allele scores' effect, --no_allele_scores beside the normal run,
was a second one 2026-10-15 to -16, until --build recorded the allele scores off.)

A row is a taxon of a sample (meta_sample, taxon_name). Each database (or run) is scored on its own rows, as the trainer
scores a table (TP: present and called; FP: absent and called; FN: present and not called; F1 = 2TP / (2TP + FP + FN)),
and with the marker truth too (machine_learning_cmdline.py's section Marker truth: a present taxon whose genomes gave no
read on a marker gene left out, the present species with such reads but no row, <table>.unseen.tsv, counted as missed).
A row's call flips when it is called in one and not in the other; a taxon with a row in one only counts as not called in
the other.

The calls of the profiled-again samples are protal's own (the dumps' prediction column: what a user's run calls) or, for
a model the run did not profile with (the other kinds of --model, or any model on a database whose own models are not the
trained ones), the model file's as protal calls with it (model_pmml.py's probabilities, the trainer's singleton veto and
the model's knob curve at the sample's depth). Where both are there, the rows whose two calls differ are counted: a parity
check on the database the users get.

pandas, numpy and the trainer are imported where they are used (build_gtdb_database.py imports this module).
"""

import collections
import os

# The counts of a side (a database or a run) on the rows compared, and of the marker truth.
SIDE_COUNTS = ("rows", "TP", "FP", "FN", "F1", "marker_TP", "marker_FN", "marker_F1")
FLIPS = ("flips_to_called_present", "flips_to_called_absent", "flips_to_not_called_present",
         "flips_to_not_called_absent")
# The noise of a refit's F1 (docs/claude: one refit moves the test F1 by ~0.001-0.002): a change beyond it is said.
F1_NOISE = 0.002


def species_name(name):
    """A species as the tables name it ("s__Genus species")."""
    name = str(name).strip()
    return name if name.startswith("s__") else "s__" + name


def read_species(path):
    """The species of a heldout_species.txt (first column), as the tables name them; an empty set without the file."""
    if not path or not os.path.isfile(path):
        return set()
    with open(path) as fh:
        return {species_name(line.split("\t")[0]) for line in fh if line.strip() and not line.startswith("#")}


def has_rows(path):
    """Whether a collector's table is there with a header at least (a collection whose samples were all left out
    writes an empty one)."""
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def _frame(path, usecols=None):
    import pandas as pd
    return pd.read_csv(path, sep="\t", low_memory=False, usecols=usecols, keep_default_na=True,
                       dtype={"meta_sample": str, "meta_scenario": str, "taxon_name": str})


def _flag(series):
    """0/1 of a truth or call column (1/0, True/False, as text or numbers)."""
    return series.astype(str).str.strip().str.lower().isin(("1", "1.0", "true")).astype(int)


SIDE_COLUMNS = ["meta_sample", "taxon_name", "truth", "call", "meta_marker_reads", "meta_novel_congener", "meta_scenario",
                "source"]
# A present species' genome in the sample (source): the database's representative, another real genome of the species
# (a real strain), or an in-silico strain (meta_rep_genome, meta_insilico_strain); "" for an absent taxon.
SOURCES = ("representative", "real strain", "in-silico strain")


def _side(frame, call):
    """The columns a side needs, its call given."""
    import numpy as np
    out = frame.reindex(columns=["meta_sample", "taxon_name", "meta_marker_reads", "meta_novel_congener",
                                 "meta_scenario"]).copy()
    out["truth"] = _flag(frame["truth"]).to_numpy()
    out["call"] = np.asarray(call, dtype=int)
    out["taxon_name"] = out["taxon_name"].map(species_name)
    out["meta_scenario"] = out["meta_scenario"].fillna("")
    present = out["truth"].to_numpy() == 1
    if "meta_rep_genome" in frame.columns:  # (a table or calls file of a collector before 2026-09-29 has none)
        rep = _flag(frame["meta_rep_genome"]).to_numpy() == 1
        insilico = (_flag(frame["meta_insilico_strain"]).to_numpy() == 1 if "meta_insilico_strain" in frame.columns
                    else np.zeros(len(frame), dtype=bool))
        source = np.where(rep, SOURCES[0], np.where(insilico, SOURCES[2], SOURCES[1]))
        out["source"] = np.where(present, source, "")
    else:
        out["source"] = ""
    return out[SIDE_COLUMNS]


def trainer_side(calls_path, samples=None):
    """The test rows of PREFIX.calls.tsv.gz (the final model's calls as protal makes them on the test samples) as a side,
    of the samples given (all without)."""
    frame = _frame(calls_path)
    frame = frame[frame["set"].astype(str) == "test"]
    if samples is not None:
        frame = frame[frame["meta_sample"].isin(samples)]
    return _side(frame, frame["call"].astype(int).to_numpy())


def model_scores(frame, model_path):
    """A model file's scores on a table's rows as protal scores and calls them by default: its probabilities
    (model_pmml.load_model, which protal's equal: check_model_parity.py) after the singleton veto (call_scores), the
    knob of its knob curve at each sample's depth (knob_at, sample_depths), or 0.5 (--knob's default) without one, and
    the calls. -> (scores, knobs, calls) arrays."""
    import numpy as np
    import machine_learning_cmdline as trainer
    from model_pmml import load_model, read_depth_knob_curve
    model = load_model(model_path)
    p = model.predict(frame)
    curve = read_depth_knob_curve(model_path)
    knobs = trainer.knob_at(curve, trainer.sample_depths(frame)) if curve else np.full(len(frame), 0.5)
    scores = np.asarray(trainer.call_scores(frame, p), dtype=float)
    return scores, np.asarray(knobs, dtype=float), (scores >= knobs).astype(int)


def model_calls(frame, model_path):
    """The calls a model file makes on a table's rows as protal makes them by default (model_scores)."""
    return model_scores(frame, model_path)[2]


# The rows whose call by protal differs from its model file's (final_check: OUT_differ.tsv beside OUT).
DIFFER_COLUMNS = ["read_type", "kind", "meta_sample", "taxon_name", "truth", "protal_probability", "protal_call",
                  "model_score", "model_knob", "model_call"]


def table_sides(table_path, models=None, protal_calls=True, differing=None):
    """A collector table's rows as sides: {"protal": protal's own calls (its prediction column)} and one per model
    file of `models` ({name: PMML path}: model_calls). -> ({name: side}, the table's samples, {name: rows whose model
    call differs from protal's} for each model, when protal_calls). differing: a dict that gets {name: those rows, as
    dicts of DIFFER_COLUMNS (but read_type), protal's probability and the model's score and knob beside the calls}."""
    frame = _frame(table_path)
    sides, differ = {}, {}
    if protal_calls and "prediction" in frame.columns:
        sides["protal"] = _side(frame, _flag(frame["prediction"]).to_numpy())
    for name, path in (models or {}).items():
        scores, knobs, calls = model_scores(frame, path)
        sides[name] = _side(frame, calls)
        if "protal" in sides:
            unlike = sides["protal"]["call"].to_numpy() != calls
            differ[name] = int(unlike.sum())
            if differing is not None:
                probability = frame["probability"].to_numpy() if "probability" in frame.columns else [None] * len(frame)
                own = sides["protal"]
                differing[name] = [{"kind": name, "meta_sample": own["meta_sample"].iat[i],
                                    "taxon_name": own["taxon_name"].iat[i], "truth": int(own["truth"].iat[i]),
                                    "protal_probability": probability[i], "protal_call": int(own["call"].iat[i]),
                                    "model_score": float(scores[i]), "model_knob": float(knobs[i]),
                                    "model_call": int(calls[i])}
                                   for i in unlike.nonzero()[0]]
    return sides, set(frame["meta_sample"].astype(str)), differ


def unseen_counts(path, samples, leave_out=frozenset()):
    """The present species without a row but with reads on a marker gene (<table>.unseen.tsv) of these samples, by
    (sample, species): a Counter; species of `leave_out` not counted."""
    out = collections.Counter()
    if not path or not os.path.isfile(path):
        return out
    frame = _frame(path)
    for sample, taxon in zip(frame["meta_sample"].astype(str), frame["taxon_name"].astype(str)):
        if sample in samples and species_name(taxon) not in leave_out:
            out[(sample, species_name(taxon))] += 1
    return out


def f1(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if tp + fp + fn else None


def side_counts(side, unseen=0):
    """SIDE_COUNTS of a side's rows; unseen: the present species without a row with marker reads (marker truth)."""
    import numpy as np
    truth, call = side["truth"].to_numpy(dtype=bool), side["call"].to_numpy(dtype=bool)
    markers = side["meta_marker_reads"].astype(float).to_numpy() if len(side) else np.zeros(0)
    blind = truth & (markers < 1)  # present without a read on a marker gene (not counted: NaN, kept)
    tp, fp, fn = int((call & truth).sum()), int((call & ~truth).sum()), int((~call & truth).sum())
    tp_m, fn_m = int((call & truth & ~blind).sum()), int((~call & truth & ~blind).sum()) + unseen
    return {"rows": len(side), "TP": tp, "FP": fp, "FN": fn, "F1": f1(tp, fp, fn), "marker_TP": tp_m,
            "marker_FN": fn_m, "marker_F1": f1(tp_m, fp, fn_m)}


def flips(before, after):
    """FLIPS: the rows (sample, taxon) of either side called in one and not in the other, by truth (a row of one side
    only is not called in the other)."""
    merged = before[["meta_sample", "taxon_name", "truth", "call"]].merge(
        after[["meta_sample", "taxon_name", "truth", "call"]], on=["meta_sample", "taxon_name"], how="outer",
        suffixes=("_before", "_after"))
    truth = merged["truth_before"].fillna(merged["truth_after"]).astype(int).to_numpy() == 1
    was = merged["call_before"].fillna(0).astype(int).to_numpy() == 1
    now = merged["call_after"].fillna(0).astype(int).to_numpy() == 1
    return {"flips_to_called_present": int((now & ~was & truth).sum()),
            "flips_to_called_absent": int((now & ~was & ~truth).sum()),
            "flips_to_not_called_present": int((was & ~now & truth).sum()),
            "flips_to_not_called_absent": int((was & ~now & ~truth).sum())}


def compare(before, after, unseen_before=None, unseen_after=None, keep=None):
    """The counts of two sides on the rows `keep` keeps (a function of a side -> a bool mask; all without), each with
    its unseen species (unseen_*: Counters by (sample, species) of those the rows kept would have, unseen_counts) for
    the marker truth: -> {before_<count>, after_<count>, F1_change (after - before), marker_F1_change, FLIPS}."""
    if keep is not None:
        before, after = before[keep(before)], after[keep(after)]
    b = side_counts(before, sum((unseen_before or {}).values()))
    a = side_counts(after, sum((unseen_after or {}).values()))
    out = {**{f"before_{k}": v for k, v in b.items()}, **{f"after_{k}": v for k, v in a.items()}}
    out["F1_change"] = difference(a["F1"], b["F1"])
    out["marker_F1_change"] = difference(a["marker_F1"], b["marker_F1"])
    out.update(flips(before, after))
    return out


def difference(a, b):
    """a - b, None if either is None; an exact 0 as 0.0 (not -0.0)."""
    if a is None or b is None:
        return None
    return a - b if a != b else 0.0


def signed(value):
    """A change with its sign, "+0.0000" for none; "-" if not known."""
    return "-" if value is None else f"{0.0 if value == 0 else value:+.4f}"


def number(value, form=".4f"):
    return "" if value is None else (format(value, form) if isinstance(value, float) else str(value))


def read_rows(path):
    """The rows of a table write_tsv wrote, as dicts: numbers as int or float, empty fields None."""
    def value(text):
        if text == "":
            return None
        for kind in (int, float):
            try:
                return kind(text)
            except ValueError:
                pass
        return text
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, (value(v) for v in line.rstrip("\n").split("\t")))) for line in fh if line.strip()]


def write_tsv(path, columns, rows):
    """rows: dicts by column; None as empty, floats to 6 significant digits."""
    with open(path + ".partial", "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join("" if row.get(c) is None else (f"{row[c]:.6g}" if isinstance(row[c], float) else str(row[c]))
                               for c in columns) + "\n")
    os.replace(path + ".partial", path)


# ---- the final check on the finished database ---------------------------------------------------------------------

FINAL_COLUMNS = (["read_type", "kind", "calls", "rows_of", "samples"] +
                 [f"{side}_{c}" for side in ("training_db", "finished_db") for c in SIDE_COUNTS] +
                 ["F1_change", "marker_F1_change"] + list(FLIPS) + ["calls_unlike_the_model_file"])
# The present rows by the genome they were simulated from (since 2026-10-16): the finished database takes its strain
# alleles from every genome of a species, the real strains the samples were simulated from among them (the training
# database leaves those out), so their gain there is partly their own alleles; the representatives' and the in-silico
# strains' is not (docs/claude/2026-10-16-r226-v31). Their FP are 0 (present rows only): read their FN and sensitivity.
SOURCE_SUBSETS = (("of them present, simulated from the representative", SOURCES[0]),
                  ("of them present, simulated from a real strain", SOURCES[1]),
                  ("of them present, simulated from an in-silico strain", SOURCES[2]))
SUBSETS = (("species not held out", None),
           ("of them without a held-out congener in the sample", "no novel congener")) + SOURCE_SUBSETS
SOURCE_RULES = {rule for _, rule in SOURCE_SUBSETS}
SOURCE_LABELS = {label for label, _ in SOURCE_SUBSETS}


def final_check(read_types, tables, unseen_tables, trainer_calls, models, default_kind, heldout, out_path):
    """model_logs/final_db_check.tsv. read_types: in order; tables: {read type: the final check's table on the finished
    database}; unseen_tables: {read type: (the test collection's <table>.unseen.tsv, the final check's)};
    trainer_calls: {(read type, kind): PREFIX.calls.tsv.gz}; models: {(read type, kind): PREFIX.xml};
    default_kind: the kind protal profiled with (the database's default models: their calls are protal's own, and
    the model file's are compared with them); heldout: the held-out species (read_species). -> (the rows written, the
    read types without samples)."""
    rows, missing = [], []
    differing_of = {}  # {read type: {kind: the rows protal calls otherwise than the model file}}
    for t in read_types:
        table = tables.get(t)
        if not has_rows(table):
            missing.append(t)
            continue
        kinds = [k for (rt, k) in trainer_calls if rt == t]
        files = {k: models[(t, k)] for k in kinds if (t, k) in models and os.path.isfile(models[(t, k)])}
        sides, samples, differ = table_sides(table, files, differing=differing_of.setdefault(t, {}))
        if not samples:
            missing.append(t)
            continue
        before_unseen_path, after_unseen_path = unseen_tables.get(t, (None, None))
        unseen_before = unseen_counts(before_unseen_path, samples, heldout)
        unseen_after = unseen_counts(after_unseen_path, samples, heldout)
        for k in kinds:
            if not os.path.isfile(trainer_calls[(t, k)]):
                continue
            after = sides.get("protal") if k == default_kind and "protal" in sides else sides.get(k)
            if after is None:
                continue
            before = trainer_side(trainer_calls[(t, k)], samples)
            for label, rule in SUBSETS:
                def keep(side, rule=rule):
                    mask = ~side["taxon_name"].isin(heldout)
                    if rule == "no novel congener":
                        mask &= side["meta_novel_congener"].astype(str).isin(("0", "0.0", "False", "false"))
                    elif rule in SOURCE_RULES:
                        mask &= side["source"] == rule
                    return mask
                if rule in SOURCE_RULES and not (after["source"] == rule).any():
                    continue  # (a table without the genome columns, or no such species)
                # the unseen species have no row: those of the subset are all the unseen (none is held out); for the
                # samples' congener and genome subsets they are left out (whether a held-out congener was there, and
                # which genome a species was simulated from, are the rows')
                counts = compare(before, after, unseen_before if rule is None else None,
                                 unseen_after if rule is None else None, keep)
                rows.append({"read_type": t, "kind": k, "calls": "protal" if after is sides.get("protal") else "model file",
                             "rows_of": label, "samples": len(samples),
                             **{f"training_db_{c}": counts[f"before_{c}"] for c in SIDE_COUNTS},
                             **{f"finished_db_{c}": counts[f"after_{c}"] for c in SIDE_COUNTS},
                             "F1_change": counts["F1_change"], "marker_F1_change": counts["marker_F1_change"],
                             **{f: counts[f] for f in FLIPS},
                             "calls_unlike_the_model_file": differ.get(k) if k == default_kind and rule is None else None})
    write_tsv(out_path, FINAL_COLUMNS, rows)
    write_tsv(differ_path(out_path), DIFFER_COLUMNS,
              [{"read_type": t, **row} for t, kinds in differing_of.items() for rows_ in kinds.values() for row in rows_])
    return rows, missing


def differ_path(out_path):
    """The table of the rows protal calls otherwise than its model file, beside final_db_check.tsv."""
    return out_path[:-len(".tsv")] + "_differ.tsv" if out_path.endswith(".tsv") else out_path + ".differ"


def final_check_lines(rows, missing, default_kind):
    """The final check in words: per read type and kind, F1 on the training database -> on the finished database over
    the rows of species not held out (and without a held-out congener), the flips, protal against the model file; and
    the warnings, those whose F1 moved beyond F1_NOISE. -> (lines, warnings)."""
    lines, warnings = [], []
    for row in rows:
        change = row["F1_change"]
        main = row["rows_of"] == SUBSETS[0][0]
        if row["rows_of"] in SOURCE_LABELS:  # present rows only: their sensitivity
            before = row["training_db_TP"] + row["training_db_FN"]
            after = row["finished_db_TP"] + row["finished_db_FN"]
            lines.append(f"{row['read_type']} {row['kind']}, {row['rows_of']}: sensitivity "
                         f"{number(row['training_db_TP'] / before if before else None)} -> "
                         f"{number(row['finished_db_TP'] / after if after else None)} (FN {row['training_db_FN']} -> "
                         f"{row['finished_db_FN']} of {before} -> {after})")
            continue
        text = (f"{row['read_type']} {row['kind']}{' (default)' if row['kind'] == default_kind else ''}, "
                f"{row['rows_of']}: F1 {number(row['training_db_F1'])} -> {number(row['finished_db_F1'])}"
                f" ({signed(change)}; FP {row['training_db_FP']} -> "
                f"{row['finished_db_FP']}, FN {row['training_db_FN']} -> {row['finished_db_FN']}" +
                (f"; marker truth F1 {number(row['training_db_marker_F1'])} -> {number(row['finished_db_marker_F1'])}"
                 if main else "") +
                f"; calls flipped: {row['flips_to_called_present'] + row['flips_to_called_absent']} to called "
                f"({row['flips_to_called_absent']} absent), {row['flips_to_not_called_present'] + row['flips_to_not_called_absent']} "
                f"to not called ({row['flips_to_not_called_present']} present))" +
                (f"; protal's calls unlike the model file's: {row['calls_unlike_the_model_file']}"
                 if row.get("calls_unlike_the_model_file") is not None else ""))
        lines.append(text)
        if main and change is not None and abs(change) > F1_NOISE:
            warnings.append(f"{row['read_type']} {row['kind']}: F1 {number(row['training_db_F1'])} on the training database, "
                            f"{number(row['finished_db_F1'])} on the finished database ({signed(change)}, beyond a "
                            f"refit's noise of {F1_NOISE}) over the {row['samples']} test samples' "
                            "rows of species not held out")
        if row.get("calls_unlike_the_model_file"):
            warnings.append(f"{row['read_type']}: protal's calls on the finished database differ from its model file's at "
                            f"{row['calls_unlike_the_model_file']} rows (a parity problem: check_model_parity.py; the rows in "
                            "final_db_check_differ.tsv)")
    if missing:
        lines.append(f"no samples written to the disk of {', '.join(missing)} (all streamed into protal, or none in the "
                     "test set): not checked")
    return lines, warnings
