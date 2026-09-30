#!/usr/bin/env python3
"""summarize.py - tables from evaluate.py outputs.

Usage: summarize.py <label> [<label> ...]   (reads results/<label>.sites.tsv etc.)
Prints, per label and stage (raw / qc), by nominal depth and control (rep genome) vs strain:
  pos       reference positions of genes present in the genome (site classes inv+snp+near+del)
  called    fraction of those holding a base or IUPAC code (not N, '-', or removed)
  acc       fraction of ACGT calls equal to the true base (inv+snp+near sites)
  iupac     fraction of calls that are IUPAC codes
  sens      true-SNP positions called with the correct alt / all true-SNP positions
  sensCov   same, among true-SNP positions that are covered (not '-' / removed)
  refAtSNP  true-SNP positions called as the reference base (per 1000 called SNP sites)
  fSNPppm   false alt calls at invariant sites per million called invariant sites
  N_inv/N_snp  N rate among covered (not '-' / removed) invariant / true-SNP positions
  errNear   ACGT error rate within 5 bp of a true indel (vs errInv at invariant sites)
"""
import os, sys
import pandas as pd

ACC = os.path.expanduser("~/audit5/accuracy")
CALLED = ["correct", "ref_at_snp", "false_alt", "wrong_alt", "iupac_with_truth", "iupac_without_truth"]
ACGTC = ["correct", "ref_at_snp", "false_alt", "wrong_alt"]


def depth_bin(d):
    for b in (1, 2, 3, 5, 10, 20, 50):
        if abs(d - b) / b < 0.1:
            return b
    return round(d, 1)


def table(df, by):
    rows = []
    for key, g in df.groupby(by):
        piv = g.pivot_table(index="site", columns="call", values="count", aggfunc="sum", fill_value=0)
        def c(site, calls):
            if site not in piv.index:
                return 0
            return int(sum(piv.loc[site].get(x, 0) for x in calls))
        sites = ["inv", "snp", "near_indel", "del"]
        tot = sum(int(piv.loc[s].sum()) for s in sites if s in piv.index)
        called = sum(c(s, CALLED) for s in ["inv", "snp", "near_indel"]) + c("del", CALLED)
        acgt = sum(c(s, ACGTC) for s in ["inv", "snp", "near_indel"])
        correct = sum(c(s, ["correct"]) for s in ["inv", "snp", "near_indel"])
        iup = sum(c(s, ["iupac_with_truth", "iupac_without_truth"]) for s in ["inv", "snp", "near_indel"])
        snp_tot = int(piv.loc["snp"].sum()) if "snp" in piv.index else 0
        snp_ok = c("snp", ["correct"])
        snp_cov = snp_tot - c("snp", ["gap", "absent"])
        inv_called = c("inv", CALLED)
        inv_cov = (int(piv.loc["inv"].sum()) if "inv" in piv.index else 0) - c("inv", ["gap", "absent"])
        near_acgt = c("near_indel", ACGTC)
        inv_acgt = c("inv", ACGTC)
        r = dict(zip(by, key if isinstance(key, tuple) else (key,)))
        r.update(pos=tot,
                 called=called / tot if tot else float("nan"),
                 acc=correct / acgt if acgt else float("nan"),
                 iupac=iup / max(called, 1),
                 sens=snp_ok / snp_tot if snp_tot else float("nan"),
                 sensCov=snp_ok / snp_cov if snp_cov else float("nan"),
                 refAtSNP=1000 * c("snp", ["ref_at_snp"]) / max(c("snp", CALLED), 1),
                 fSNPppm=1e6 * c("inv", ["false_alt"]) / max(inv_called, 1),
                 N_inv=c("inv", ["N"]) / inv_cov if inv_cov else float("nan"),
                 N_snp=c("snp", ["N"]) / snp_cov if snp_cov else float("nan"),
                 errInv=1 - c("inv", ["correct"]) / inv_acgt if inv_acgt else float("nan"),
                 errNear=1 - c("near_indel", ["correct"]) / near_acgt if near_acgt else float("nan"),
                 n_snp=snp_tot)
        rows.append(r)
    return pd.DataFrame(rows)


def load(label):
    df = pd.read_csv(os.path.join(ACC, "results", label + ".sites.tsv"), sep="\t")
    df = df[df.site != "absent_gene"].copy()
    df["depth"] = df.depth.map(depth_bin)
    df["kind"] = df.genome.map(lambda g: "rep" if g.startswith("GCF") else "strain")
    df["domain"] = df.species.map(lambda s: "archaea" if s.startswith("Calidella") else "bacteria")
    return df


def fmt(t):
    f = {"called": "{:.3f}", "acc": "{:.5f}", "iupac": "{:.4f}", "sens": "{:.3f}", "sensCov": "{:.3f}",
         "refAtSNP": "{:.1f}", "fSNPppm": "{:.0f}", "N_inv": "{:.4f}", "N_snp": "{:.4f}", "errInv": "{:.5f}",
         "errNear": "{:.4f}"}
    t = t.copy()
    for k, v in f.items():
        if k in t:
            t[k] = t[k].map(lambda x: v.format(x) if pd.notna(x) else "nan")
    return t.to_string(index=False)


def main():
    pd.set_option("display.width", 250)
    for label in sys.argv[1:]:
        df = load(label)
        print(f"=== {label}")
        for stage in ("raw", "qc"):
            d = df[df.stage == stage]
            print(f"--- stage {stage}: by depth x kind")
            print(fmt(table(d, ["depth", "kind"])))
        print("--- raw by domain (strains)")
        print(fmt(table(df[(df.stage == "raw") & (df.kind == "strain")], ["domain", "depth"])))
        rows = pd.read_csv(os.path.join(ACC, "results", label + ".rows.tsv"), sep="\t")
        rows["depth"] = rows.depth.map(lambda d: depth_bin(float(d)))
        print("--- row status (raw, qc) by depth")
        print(rows.groupby(["depth", "raw_status", "qc_status"]).size().unstack(fill_value=0).to_string())


def indel_tables(label):
    """Per true indel event: status by depth (raw / qc); false indels per sample row."""
    ind = pd.read_csv(os.path.join(ACC, "results", label + ".indels.tsv"), sep="\t")
    ind["depth"] = ind.depth.map(depth_bin)
    ind = ind[~ind.genome.str.startswith("GCF")]
    for stage in ("raw", "qc"):
        d = ind[ind.stage == stage]
        t = d.groupby(["type", "depth", "status"]).size().unstack(fill_value=0)
        t["n"] = t.sum(axis=1)
        print(f"--- {label} true indel events, stage {stage}")
        print(t.to_string())
    fi = pd.read_csv(os.path.join(ACC, "results", label + ".falseindels.tsv"), sep="\t")
    fi["depth"] = fi.depth.map(depth_bin)
    print(f"--- {label} false indels (raw): insertion columns filled where truth has none; "
          f"internal '-' runs flanked by calls where truth has no deletion")
    print(fi[fi.stage == "raw"].groupby(["type", "depth"]).agg(n=("pos", "size"), mean_len=("len", "mean"),
                                                               len3=("len", lambda x: (x == 3).sum())).to_string())


def mix_tables(label):
    """Mixtures: what the row holds where the two strains differ, by minor fraction and total depth."""
    m = pd.read_csv(os.path.join(ACC, "results", label + ".mix.tsv"), sep="\t")
    m["depth"] = m.depth.round()
    m["diff"] = m.site.isin(["major_alt_only", "minor_alt_only", "both_alt_differ"])
    order = ["major_base", "minor_base", "iupac_both", "iupac_other", "N", "gap", "absent", "other"]
    for stage in ("raw", "qc"):
        d = m[(m.stage == stage) & m["diff"]]
        t = d.pivot_table(index=["depth", "minor_frac"], columns="call", values="count", aggfunc="sum", fill_value=0)
        t = t.reindex(columns=[c for c in order if c in t.columns])
        n = t.sum(axis=1)
        f = t.div(n, axis=0).round(3)
        f["n_sites"] = n
        print(f"--- {label} stage {stage}: calls at sites where the two strains differ (fraction)")
        print(f.to_string())
        # among positions holding a base/IUPAC: fraction reporting the mixture
        cov = t.drop(columns=[c for c in ("N", "gap", "absent") if c in t.columns]).sum(axis=1)
        print("  of called:", (t.get("iupac_both", 0) / cov).round(3).to_dict())
    b = m[(m.stage == "raw") & (m.site == "both_ref")]
    t = b.pivot_table(index=["depth", "minor_frac"], columns="call", values="count", aggfunc="sum", fill_value=0)
    print(f"--- {label} raw: sites where both strains carry the reference base")
    print(t.to_string())


if __name__ == "__main__":
    if os.environ.get("MIX"):
        for label in sys.argv[1:]:
            mix_tables(label)
    elif os.environ.get("INDELS"):
        for label in sys.argv[1:]:
            indel_tables(label)
    else:
        main()
