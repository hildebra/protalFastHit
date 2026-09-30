#!/usr/bin/env python3
"""evaluate.py - (copied from ../../2026-09-29-strain-audit/scripts/accuracy/evaluate.py; the only change: the
instrumented qcmsa is $QCMSA_COLMAP, so that each version's qcmsa is re-run) score protal's strain MSAs against the true genome sequences.

Usage: evaluate.py <run> <protal_out_dir> <label>
  <run>   design manifest name (A, B or C) giving each sample's genomes and depths
  <label> free text for the output files (results/<label>.*.tsv)

For every species MSA in <out>/strains: map raw-MSA columns to (gene, reference position) via the
partition file and the reference row ('-' in the reference row = insertion column placed before
the next reference base), check the reference row against the representative's gene, and classify
every (sample row, reference position) by the true base of that sample's genome (truth.pkl) and
by what the row holds. The qcmsa-filtered MSA is mapped back with qcmsa_colmap.py (verbatim
qcmsa.py + column map), whose output is checked to be identical to protal's <species>.msa.fna.

Writes results/<label>.sites.tsv (counts by site class x call class per sample row and stage),
results/<label>.indels.tsv (per true indel event), results/<label>.falseindels.tsv,
results/<label>.rows.tsv (row status per sample x species) and, for mixtures (run B),
results/<label>.mix.tsv.
"""
import os, pickle, shlex, subprocess, sys, tempfile, re
from collections import defaultdict
import numpy as np

ACC = os.path.expanduser("~/audit5/accuracy")
TRUTH = pickle.load(open(os.path.join(ACC, "truth.pkl"), "rb"))
GENES, TR, META = TRUTH["genes"], TRUTH["truth"], TRUTH["meta"]
REP_OF = META["rep_of"]
WIN = 5  # bases either side of a true indel counted as 'near_indel'

IUPAC = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
B = lambda c: ord(c)
GAP, NN, XX = B("-"), B("N"), B("X")
ACGT = np.zeros(256, bool)
for c in "ACGT":
    ACGT[B(c)] = True
IUP = np.zeros(256, bool)
for c in IUPAC:
    IUP[B(c)] = True
# IUPAC membership table: IN[code, base] = base is in code (for ACGT codes: identity)
IN = np.zeros((256, 256), bool)
for c in "ACGT":
    IN[B(c), B(c)] = True
for c, s in IUPAC.items():
    for b in s:
        IN[B(c), B(b)] = True

SITE = ["inv", "snp", "del", "near_indel"]
CALL = ["correct", "ref_at_snp", "false_alt", "wrong_alt", "iupac_with_truth", "iupac_without_truth", "N", "gap", "absent"]


def read_fasta(path):
    names, seqs, cur = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if names:
                    seqs.append("".join(cur))
                names.append(line[1:].strip())
                cur = []
            else:
                cur.append(line)
    if names:
        seqs.append("".join(cur))
    return names, seqs


# [v07-vs-v06] 0.6.0a writes 0-based partitions (fixed in round 2): PARTITION_OFFSET=1 reads them.
PART_OFFSET = int(os.environ.get("PARTITION_OFFSET", "0"))


def read_partition(path):
    parts = []
    for line in open(path):
        m = re.search(r"gene(\d+)\s*=\s*(\d+)-(\d+)", line)
        if m:
            parts.append((int(m.group(1)), int(m.group(2)) - 1 + PART_OFFSET, int(m.group(3)) - 1 + PART_OFFSET))
    return parts


def load_manifest(run):
    d = defaultdict(lambda: defaultdict(list))
    with open(os.path.join(ACC, "design", run + ".manifest.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            d[f["sample"]][f["species"]].append((f["genome"], float(f["vertical_coverage"])))
    return d


def qc_rows(strains, name, log_path, tmpdir):
    """Re-run qcmsa (instrumented copy) exactly as protal did; return {sample: full-length row
    array with 'X' in removed columns} or None if protal wrote no filtered MSA."""
    final = os.path.join(strains, name + ".msa.fna")
    if not os.path.exists(final):
        return None, "no_qcmsa_output"
    cmd = None
    raw = os.path.join(strains, name + ".raw.msa.fna")
    for line in open(log_path, errors="replace"):
        i = line.find("[qcmsa] ")
        if i >= 0 and ("'" + raw + "'") in line:
            cmd = shlex.split(line[i + 8:].strip())
    if cmd is None:
        raise SystemExit(f"no qcmsa command for {name} in {log_path}")
    # cmd: [python3] script msa part meta --prefix P --reapply-hcov N [extra...]
    if cmd[0] == "python3":
        cmd = cmd[1:]
    args = cmd[1:]
    k = args.index("--prefix")
    prefix = os.path.join(tmpdir, name)
    args[k + 1] = prefix
    # [v07-vs-v06] the instrumented copy of the qcmsa the run used (make_colmap.py)
    colmap = os.environ.get("QCMSA_COLMAP", os.path.join(ACC, "qcmsa_colmap.py"))
    subprocess.run(["python3", colmap] + args, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if open(prefix + ".msa.fna").read() != open(final).read():
        raise SystemExit(f"instrumented qcmsa output differs from protal's for {name}")
    with open(prefix + ".colmap.txt") as fh:
        kept = fh.readline().rstrip("\n").split("\t")[1:]
        cols = np.array([int(x) for x in fh.read().split()], dtype=np.int64)
    names, seqs = read_fasta(final)
    assert names == kept
    return dict(names=names, seqs=seqs, cols=cols), "ok"


def classify(ref, tb, calls):
    """ref, tb, calls: uint8 arrays over reference positions. Returns (site, call) code arrays."""
    L = len(ref)
    site = np.zeros(L, np.int8)
    snp = (tb != ref) & (tb != GAP)
    site[snp] = 1
    site[tb == GAP] = 2
    call = np.full(L, -1, np.int8)
    acgt = ACGT[calls]
    call[acgt & (calls == tb)] = 0
    call[acgt & (calls != tb) & (calls == ref) & snp] = 1
    call[acgt & (calls != tb) & (calls != ref) & (tb == ref)] = 2
    call[acgt & (calls != tb) & (call < 0)] = 3
    iup = IUP[calls]
    call[iup & IN[calls, tb]] = 4
    call[iup & ~IN[calls, tb]] = 5
    call[calls == NN] = 6
    call[calls == GAP] = 7
    call[calls == XX] = 8
    assert (call >= 0).all(), set(calls[call < 0].tolist())
    return site, call


def near_mask(L, t):
    m = np.zeros(L, bool)
    tb = t["tb"]
    p = 0
    while p < L:
        if tb[p] == "-":
            q = p
            while q < L and tb[q] == "-":
                q += 1
            m[max(0, p - WIN):min(L, q + WIN)] = True
            p = q
        else:
            p += 1
    for p in t["ins"]:
        m[max(0, p - WIN):min(L, p + WIN)] = True
    return m


def main():
    run, out, label = sys.argv[1], os.path.abspath(sys.argv[2]), sys.argv[3]
    log_path = sys.argv[4] if len(sys.argv) > 4 else None
    man = load_manifest(run)
    strains = os.path.join(out, "strains")
    os.makedirs(os.path.join(ACC, "results"), exist_ok=True)
    tmpdir = tempfile.mkdtemp(prefix="qc_", dir=os.path.join(ACC, "results"))
    site_counts = defaultdict(int)
    indel_rows, false_indels, row_rows, mix_counts = [], [], [], defaultdict(int)
    error_rows = []
    refcheck = defaultdict(int)
    species_all = sorted({sp for s in man.values() for sp in s})
    for sp in species_all:
        name = "s__" + sp.replace(" ", "_")
        rawp = os.path.join(strains, name + ".raw.msa.fna")
        has_msa = os.path.exists(rawp)
        if has_msa:
            names, seqs = read_fasta(rawp)
            parts = read_partition(os.path.join(strains, name + ".raw.partition.txt"))
            assert names[0] == name + "_reference"
            refrow = np.frombuffer(seqs[0].encode(), np.uint8)
            rows = {n: np.frombuffer(s.encode(), np.uint8) for n, s in zip(names[1:], seqs[1:])}
            qc, qc_status = qc_rows(strains, name, log_path, tmpdir)
            qcrows = {}
            if qc is not None:
                for n, s in zip(qc["names"], qc["seqs"]):
                    full = np.full(len(refrow), XX, np.uint8)
                    full[qc["cols"]] = np.frombuffer(s.encode(), np.uint8)
                    qcrows[n] = full
            meta_samples = set()
            for line in open(os.path.join(strains, name + ".meta.tsv")):
                meta_samples.add(line.split("\t")[0])
        else:
            names, parts, rows, qcrows, meta_samples, qc_status = [name + "_reference"], [], {}, {}, set(), "no_msa"
        # gene -> (refcols, inscols, ins_before_pos)
        gmap = {}
        for gid, s, e in parts:
            cols = np.arange(s, e + 1)
            seg = refrow[cols]
            refcols = cols[seg != GAP]
            inscols = cols[seg == GAP]
            # reference position each insertion column precedes
            ins_pos = np.searchsorted(refcols, inscols)
            gref = GENES.get((sp, gid))
            ok = gref is not None and bytes(refrow[refcols]).decode() == gref
            refcheck["ok" if ok else "bad"] += 1
            gmap[gid] = (refcols, inscols, ins_pos)
        species_genes = sorted(g for (s2, g) in GENES if s2 == sp)
        for sample in sorted(man):
            gl = man[sample].get(sp)
            if not gl:
                continue
            present = sample in rows
            status = "present" if present else ("hcov_dropped" if sample in meta_samples else "not_profiled")
            qc_present = sample in qcrows
            depth = sum(d for _, d in gl)
            row_rows.append((sample, sp, ";".join(g for g, _ in gl), "%.3g" % depth, status,
                             "present" if qc_present else ("removed" if present and qc_status == "ok" else qc_status)))
            for stage in ("raw", "qc"):
                row = rows.get(sample) if stage == "raw" else qcrows.get(sample)
                for gid in species_genes:
                    ref = np.frombuffer(GENES[(sp, gid)].encode(), np.uint8)
                    L = len(ref)
                    if row is not None and gid in gmap:
                        refcols, inscols, ins_pos = gmap[gid]
                        calls = row[refcols]
                    else:
                        calls = np.full(L, XX, np.uint8)
                        inscols, ins_pos = np.array([], np.int64), np.array([], np.int64)
                    if len(gl) == 1:
                        genome = gl[0][0]
                        t = TR.get((genome, gid))
                        if t is None:
                            # gene absent from this genome: anything called here came from elsewhere
                            nc = int(np.sum(ACGT[calls] | IUP[calls]))
                            site_counts[(stage, sample, sp, genome, gl[0][1], "absent_gene", "called")] += nc
                            site_counts[(stage, sample, sp, genome, gl[0][1], "absent_gene", "not_called")] += L - nc
                            continue
                        tb = np.frombuffer("".join(t["tb"]).encode(), np.uint8)
                        site, call = classify(ref, tb, calls)
                        near = near_mask(L, t)
                        site[near & (site != 2)] = 3
                        key = site.astype(np.int64) * len(CALL) + call
                        bc = np.bincount(key, minlength=len(SITE) * len(CALL))
                        for idx in np.nonzero(bc)[0]:
                            site_counts[(stage, sample, sp, genome, gl[0][1], SITE[idx // len(CALL)], CALL[idx % len(CALL)])] += int(bc[idx])
                        # --- positions worth a pileup: wrong/ambiguous calls, and N at true SNPs
                        if stage == "raw":
                            sel = np.isin(call, [1, 2, 3, 4, 5]) | ((call == 6) & (site == 1))
                            for p in np.nonzero(sel & (site != 2))[0]:
                                error_rows.append((sample, sp, genome, "%.3g" % gl[0][1], gid, int(p), SITE[site[p]],
                                                   CALL[call[p]], chr(calls[p]), chr(tb[p]), chr(ref[p])))
                        # --- indel events
                        called = ACGT[calls] | IUP[calls]
                        rowins = defaultdict(list)
                        if row is not None and len(inscols):
                            ich = row[inscols]
                            for p, c in zip(ins_pos, ich):
                                if c != GAP and c != XX:
                                    rowins[int(p)].append(chr(c))
                        tbl = t["tb"]
                        p = 0
                        true_del_starts = set()
                        while p < L:
                            if tbl[p] == "-":
                                q = p
                                while q < L and tbl[q] == "-":
                                    q += 1
                                true_del_starts.update(range(max(0, p - WIN), min(L, q + WIN)))
                                cells = calls[p:q]
                                fl = (p == 0 or called[p - 1]) and (q >= L or called[q])
                                if (cells == XX).all():
                                    st = "absent"
                                elif (cells == GAP).all():
                                    st = "represented" if fl else "gap_no_flank"
                                elif ACGT[cells].any():
                                    st = "missed_bases_called"
                                elif (cells == NN).any():
                                    st = "N"
                                else:
                                    st = "other"
                                if st == "missed_bases_called":
                                    # shifted representation: an equal-length gap run nearby
                                    lo, hi = max(0, p - WIN), min(L, q + WIN)
                                    w = calls[lo:hi]
                                    s_ = "".join("-" if c == GAP else "." for c in w)
                                    if ("." + "-" * (q - p) + ".") in ("." + s_ + "."):
                                        st = "represented_shifted"
                                indel_rows.append((stage, sample, sp, genome, "%.3g" % gl[0][1], gid, "del", p, q - p, st))
                                p = q
                            else:
                                p += 1
                        for p, s in t["ins"].items():
                            got = "".join(rowins.get(p, []))
                            fl = (p == 0 or called[p - 1]) and (p >= L or called[p])
                            if (calls[max(0, p - 1):p + 1] == XX).all():
                                st = "absent"
                            elif got == s:
                                st = "represented"
                            elif got:
                                st = "wrong_sequence"
                            else:
                                near_ins = [q for q in rowins if abs(q - p) <= WIN and len("".join(rowins[q])) == len(s)]
                                st = "represented_shifted" if near_ins else ("missed" if fl else "no_flank_cov")
                            indel_rows.append((stage, sample, sp, genome, "%.3g" % gl[0][1], gid, "ins", p, len(s), st))
                        # --- false indels in the row
                        true_ins_near = set()
                        for p in t["ins"]:
                            true_ins_near.update(range(p - WIN, p + WIN + 1))
                        for p, chars in rowins.items():
                            if p not in true_ins_near:
                                false_indels.append((stage, sample, sp, genome, "%.3g" % gl[0][1], gid, "ins", p, len(chars)))
                        # internal '-' runs flanked by called bases (deletion or coverage dip: indistinguishable)
                        g = calls == GAP
                        if g.any():
                            idx = np.nonzero(g)[0]
                            starts = idx[np.r_[True, np.diff(idx) > 1]]
                            ends = idx[np.r_[np.diff(idx) > 1, True]] + 1
                            for a, b in zip(starts, ends):
                                if a > 0 and b < L and called[a - 1] and called[b] and a not in true_del_starts:
                                    false_indels.append((stage, sample, sp, genome, "%.3g" % gl[0][1], gid, "gap_run", int(a), int(b - a)))
                    else:
                        # mixture: two genomes (major first)
                        (g1, d1), (g2, d2) = gl[0], gl[1]
                        t1, t2 = TR.get((g1, gid)), TR.get((g2, gid))
                        if t1 is None or t2 is None:
                            continue
                        tb1 = np.frombuffer("".join(t1["tb"]).encode(), np.uint8)
                        tb2 = np.frombuffer("".join(t2["tb"]).encode(), np.uint8)
                        near = near_mask(L, t1) | near_mask(L, t2) | (tb1 == GAP) | (tb2 == GAP)
                        cls = np.full(L, "", object)
                        both_ref = (tb1 == ref) & (tb2 == ref)
                        shared = (tb1 == tb2) & (tb1 != ref)
                        maj = (tb1 != ref) & (tb2 == ref)
                        mino = (tb2 != ref) & (tb1 == ref)
                        diff2 = (tb1 != tb2) & (tb1 != ref) & (tb2 != ref)
                        fmin = d2 / (d1 + d2)
                        for cname, mask in (("both_ref", both_ref), ("shared_alt", shared), ("major_alt_only", maj),
                                            ("minor_alt_only", mino), ("both_alt_differ", diff2)):
                            mask = mask & ~near
                            if not mask.any():
                                continue
                            c = calls[mask]
                            a1, a2 = tb1[mask], tb2[mask]
                            acgt = ACGT[c]
                            iup = IUP[c]
                            res = np.full(c.shape, "other", object)
                            res[acgt & (c == a1)] = "major_base"
                            res[acgt & (c == a2) & (a1 != a2)] = "minor_base"
                            res[iup & IN[c, a1] & IN[c, a2] & (a1 != a2)] = "iupac_both"
                            res[iup & ~(IN[c, a1] & IN[c, a2])] = "iupac_other"
                            res[c == NN] = "N"
                            res[c == GAP] = "gap"
                            res[c == XX] = "absent"
                            for r_, n_ in zip(*np.unique(res, return_counts=True)):
                                mix_counts[(stage, sample, sp, g1, g2, "%.3g" % (d1 + d2), "%.3g" % fmin, cname, r_)] += int(n_)
    lab = os.path.join(ACC, "results", label)
    with open(lab + ".sites.tsv", "w") as fh:
        fh.write("stage\tsample\tspecies\tgenome\tdepth\tsite\tcall\tcount\n")
        for k, v in sorted(site_counts.items()):
            fh.write("\t".join(str(x) for x in k) + "\t%d\n" % v)
    with open(lab + ".indels.tsv", "w") as fh:
        fh.write("stage\tsample\tspecies\tgenome\tdepth\tgene\ttype\tpos\tlen\tstatus\n")
        for r in indel_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(lab + ".falseindels.tsv", "w") as fh:
        fh.write("stage\tsample\tspecies\tgenome\tdepth\tgene\ttype\tpos\tlen\n")
        for r in false_indels:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(lab + ".errors.tsv", "w") as fh:
        fh.write("sample\tspecies\tgenome\tdepth\tgene\tpos\tsite\tcall\tcalled\ttrue\tref\n")
        for r in error_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(lab + ".rows.tsv", "w") as fh:
        fh.write("sample\tspecies\tgenomes\tdepth\traw_status\tqc_status\n")
        for r in row_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    if mix_counts:
        with open(lab + ".mix.tsv", "w") as fh:
            fh.write("stage\tsample\tspecies\tmajor\tminor\tdepth\tminor_frac\tsite\tcall\tcount\n")
            for k, v in sorted(mix_counts.items()):
                fh.write("\t".join(str(x) for x in k) + "\t%d\n" % v)
    print(f"{label}: reference-row check per gene: {dict(refcheck)}", file=sys.stderr)


if __name__ == "__main__":
    main()
