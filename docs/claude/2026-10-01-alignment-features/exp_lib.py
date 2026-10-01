"""Shared code of the alignment-feature experiments on the V2 build: tables, SAM-derived features, the forest
(as scripts/random_forest_cmdline.py fits it), species-held-out cross-validation and the test set."""
import collections
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "scripts"))
from model_features import NORMALIZED_FEATURES  # noqa: E402

V2 = os.path.expanduser(os.environ.get("V2", "~/tune/V2"))
EXP = os.path.expanduser(os.environ.get("FPEXP", "~/fpexp"))
READ_TYPES = ("pe", "se", "pb", "ont")
TABLE = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv", "ont": "training_data_ont.tsv"}
KNOB = 0.5


# ---- tables -------------------------------------------------------------------------------------------------

def load(split, rt):
    """The training or test table of a read type, with genus and sample identity."""
    df = pd.read_csv(f"{V2}/{split}/{TABLE[rt]}", sep="\t", float_precision="round_trip", low_memory=False)
    df["truth"] = df["truth"].astype(int)
    df["genus"] = df["taxon_name"].str.split(" ").str[0].str[3:]
    return df


def sample_sams(split):
    """SAMPLEID -> SAM path of the run's map (profile_all/samples.map)."""
    out = {}
    with open(f"{V2}/{split}/profile_all/samples.map") as fh:
        header = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "#SAMPLEID":
                header = [x.lstrip("#") for x in f]
            elif not f[0].startswith("#") and header:
                r = dict(zip(header, f))
                out[r["SAMPLEID"]] = r["SAM"]
    return out


def simulated_genera(split):
    """sample -> set of genera simulated in it (present or held out), for scoring genus-level reports."""
    by_comm = {}
    for design in os.listdir(f"{V2}/{split}/points"):
        path = f"{V2}/{split}/points/{design}/sim/manifest.tsv"
        if os.path.exists(path):
            m = pd.read_csv(path, sep="\t")
            for sample, g in m.groupby("sample"):
                by_comm[sample] = {t.split(";")[5][3:] for t in g["taxonomy"]}
    out = dict(by_comm)
    for sample in list(by_comm):
        out[sample + "_se"] = by_comm[sample]
    for design in os.listdir(f"{V2}/{split}/points"):
        path = f"{V2}/{split}/points/{design}/sim/samples.tsv"
        if os.path.exists(path):
            for _, r in pd.read_csv(path, sep="\t").iterrows():
                out[r["sample"]] = by_comm[r["community"]]
    return out


# ---- SAM-derived features -----------------------------------------------------------------------------------

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_counts(cigar):
    c = collections.Counter()
    for n, op in CIGAR.findall(cigar):
        c[op] += int(n)
    return c


def sam_records(path):
    cat = "zstdcat" if path.endswith(".zst") else "zcat" if path.endswith(".gz") else "cat"
    proc = subprocess.Popen(f"{cat} '{path}' | grep -v '^@'", shell=True, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        f = line.split("\t", 11)
        yield f
    proc.wait()


def primary_features(split, rt, cache=True):
    """Per (sample, taxid) from the run's own SAMs (one alignment per read): mate concordance, MAPQ and clipping."""
    out_path = f"{EXP}/features/{split}_{rt}_primary.tsv"
    if cache and os.path.exists(out_path):
        return pd.read_csv(out_path, sep="\t")
    rows = []
    tables = load(split, rt)
    sams = sample_sams(split)
    for sample in sorted(tables["meta_sample"].unique()):
        acc = collections.defaultdict(lambda: collections.Counter())
        frag_both = collections.defaultdict(set)
        frag_all = collections.defaultdict(set)
        for f in sam_records(sams[sample]):
            flag = int(f[1])
            if flag & 4 or flag & 256:
                continue
            taxid = f[2].split("_")[0]
            mapq, cig = int(f[4]), cigar_counts(f[5])
            a = acc[taxid]
            a["records"] += 1
            a["low_mapq"] += mapq < 10
            a["mapq0"] += mapq == 0
            clip = cig["S"] + cig["H"]
            length = cig["M"] + cig["="] + cig["X"] + cig["I"] + clip
            a["clipped"] += clip >= max(10, 0.1 * length)
            a["clip_bases"] += clip
            a["read_bases"] += length
            frag_all[taxid].add(f[0])
            if flag & 1 and not flag & 8 and (f[6] == "=" or f[6].split("_")[0] == taxid):
                frag_both[taxid].add(f[0])
        for taxid, a in acc.items():
            n = max(1, a["records"])
            rows.append({"meta_sample": sample, "taxon": int(taxid), "sam_records": a["records"],
                         "low_mapq_share": a["low_mapq"] / n, "mapq0_share": a["mapq0"] / n,
                         "clipped_share": a["clipped"] / n, "clip_base_share": a["clip_bases"] / max(1, a["read_bases"]),
                         "mate_concordance": (len(frag_both[taxid]) / len(frag_all[taxid])) if rt == "pe" else 1.0})
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    return df


def db_genera():
    """taxid (as text) -> genus of the training database's species, from its genome2tiid.tsv."""
    out = {}
    with open(f"{V2}/training_db/genome2tiid.tsv") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                out[f[1]] = f[3].split(";")[5][3:]
    return out


def secondary_features(split, rt, deltas=(0, 1, 2), cache=True):
    """Per (sample, taxid) from the realignment with -m 3 (~/fpexp/<split>/sam): the share of the taxon's reads
    whose alternative alignment to another species of the same genus (or of another genus) is within d edits of the
    best one. Edits: mismatches + inserted + deleted bases + unaligned (clipped) read bases, summed over the mates."""
    out_path = f"{EXP}/features/{split}_{rt}_secondary.tsv"
    if cache and os.path.exists(out_path):
        return pd.read_csv(out_path, sep="\t")
    taxon_genus = db_genera()
    names = sample_sams(split)
    tables = load(split, rt)
    rows = []
    for sample in sorted(tables["meta_sample"].unique()):
        path = f"{EXP}/{split}/sam/{os.path.basename(names[sample])}"
        # Each primary or supplementary record (a mate, a single read, a long read's segment) with the secondary
        # records that follow it for the same read and mate: [taxid, edits, [(taxid, edits), ...]].
        units, last = [], {}
        for f in sam_records(path):
            flag = int(f[1])
            if flag & 4:
                continue
            c = cigar_counts(f[5])
            rec = (f[2].split("_")[0], c["X"] + c["I"] + c["D"] + c["S"] + c["H"])
            key = (f[0], flag & 0xC0)
            if flag & 256:
                if key in last:
                    last[key][2].append(rec)
            else:
                unit = [rec[0], rec[1], []]
                units.append(unit)
                last[key] = unit
        per = collections.defaultdict(collections.Counter)
        for taxid, edits, alts in units:
            g = taxon_genus.get(taxid)
            c = per[taxid]
            c["reads"] += 1
            for delta in deltas:
                c[f"congener_within_{delta}"] += any(t != taxid and taxon_genus.get(t) == g and e <= edits + delta for t, e in alts)
                c[f"other_genus_within_{delta}"] += any(taxon_genus.get(t) != g and e <= edits + delta for t, e in alts)
        for taxid, c in per.items():
            r = {"meta_sample": sample, "taxon": int(taxid), "realigned_reads": c["reads"]}
            for delta in deltas:
                r[f"congener_share_d{delta}"] = c[f"congener_within_{delta}"] / c["reads"]
                r[f"other_genus_share_d{delta}"] = c[f"other_genus_within_{delta}"] / c["reads"]
            rows.append(r)
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    return df


# ---- the depth-scaled identity ------------------------------------------------------------------------------

def identity_model(train):
    """mu, sigma_between, sigma_read of present taxa's identity: var(identity | n fragments) = sb^2 + sr^2 / n,
    fitted on the training table's present taxa by bins of n (method of moments)."""
    tp = train[(train["truth"] == 1) & (train["fragments"] > 0)].copy()
    mu = float(tp.loc[tp["fragments"] >= 20, "identity"].median())
    tp["bin"] = pd.qcut(np.log(tp["fragments"]), 10, duplicates="drop")
    g = tp.groupby("bin", observed=True).agg(v=("identity", "var"), n=("fragments", "median")).dropna()
    A = np.vstack([np.ones(len(g)), 1 / g["n"].to_numpy()]).T
    sb2, sr2 = np.linalg.lstsq(A, g["v"].to_numpy(), rcond=None)[0]
    return mu, float(np.sqrt(max(sb2, 1e-8))), float(np.sqrt(max(sr2, 1e-8)))


def identity_z(df, model):
    mu, sb, sr = model
    n = df["fragments"].clip(lower=1)
    return (df["identity"] - mu) / np.sqrt(sb ** 2 + sr ** 2 / n)


# ---- the forest, folds and scores ---------------------------------------------------------------------------

def forest(seed, n_jobs=6):
    return RandomForestClassifier(n_estimators=64, max_leaf_nodes=128, min_samples_leaf=1, max_features="sqrt",
                                  class_weight="balanced", random_state=seed, n_jobs=n_jobs)


def cv_predict(train, features, seed):
    """Species held out: GroupKFold(5) by taxon, as the trainer's 'species' scheme."""
    X, y = train[features].to_numpy(float), train["truth"].to_numpy()
    p = np.full(len(y), np.nan)
    for tr, te in GroupKFold(5).split(X, y, train["taxon"].astype(str)):
        p[te] = forest(seed).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return p


def test_predict(train, test, features, seed):
    rf = forest(seed).fit(train[features].to_numpy(float), train["truth"].to_numpy())
    return rf.predict_proba(test[features].to_numpy(float))[:, 1]


def scores(y, call):
    tp, fp, fn = int((call & (y == 1)).sum()), int((call & (y == 0)).sum()), int((~call & (y == 1)).sum())
    return {"TP": tp, "FP": fp, "FN": fn, "F1": 2 * tp / max(1, 2 * tp + fp + fn)}


def bootstrap_diff(df, call_a, call_b, n=2000, seed=1):
    """F1(b) - F1(a) on the test set with samples resampled: mean and 95% interval."""
    rng = np.random.default_rng(seed)
    samples = df["meta_sample"].to_numpy()
    uniq = np.unique(samples)
    y = df["truth"].to_numpy()
    idx = {s: np.flatnonzero(samples == s) for s in uniq}
    diffs = []
    for _ in range(n):
        pick = np.concatenate([idx[s] for s in rng.choice(uniq, len(uniq))])
        diffs.append(scores(y[pick], call_b[pick])["F1"] - scores(y[pick], call_a[pick])["F1"])
    d = np.array(diffs)
    return float(d.mean()), float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))
