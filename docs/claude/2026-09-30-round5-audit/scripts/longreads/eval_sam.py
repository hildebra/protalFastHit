"""Checks protal's long-read SAM against the truth of gen_reads.py.

  eval_sam.py SAM READS.fq TRUTH.pkl --gtdb G --db D [--fasta] [--show N]

Per record: hard clips + SEQ = read length, SEQ = the read's bases (reverse-complemented for 0x10),
QUAL = the read's qualities (reversed for 0x10), the taxon is the read's species, and the CIGAR's
aligned bases sit where the read's bases truly come from (exact on reference genomes, within 30 bp
on strains). Per read: one primary; every marker gene wholly on the read found once (primary or
supplementary), partial ones counted; records of genes not on the read."""
import argparse
import collections
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrlib import World, parse_sam, cigar_ops, revcomp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sam")
    ap.add_argument("reads")
    ap.add_argument("truth")
    ap.add_argument("--gtdb", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--fasta_qual", default=None, help="expected constant quality char (FASTA input)")
    ap.add_argument("--show", type=int, default=5)
    ap.add_argument("--full", type=float, default=0.99)
    ap.add_argument("--exact", type=float, default=0.98)
    a = ap.parse_args()
    world = World(a.gtdb, a.db)
    T = pickle.load(open(a.truth, "rb"))
    truth = {t["name"]: t for t in T["reads"]}
    offsets = T["contig_offsets"]
    seqs = {}
    with open(a.reads) as fh:
        lines = fh.read().split("\n")
    if lines[0].startswith("@"):
        for i in range(0, len(lines) - 3, 4):
            seqs[lines[i][1:].split()[0]] = (lines[i + 1], lines[i + 3])
    else:
        name = None
        for l in lines:
            if l.startswith(">"):
                name = l[1:].split()[0]
            elif name:
                seqs[name] = (l, None)
    # genes per genome in chromosome coordinates: geneid -> (start, end, strand)
    genes = {}
    for acc in offsets:
        g = {}
        for marker, contig, s, e, strand in world.markers.get(acc, []):
            gid = world.geneid[marker]
            g.setdefault(gid, []).append((offsets[acc][contig] + s, offsets[acc][contig] + e, strand))
        genes[acc] = g

    cache = {"name": None}

    def origin_of(name):
        if cache["name"] != name:
            t = truth[name]
            o = [-1] * t["read_length"]
            step = 1 if t["strand"] == "+" else -1
            for rs, cs, ln in t["blocks"]:
                for k in range(ln):
                    o[rs + k] = cs + step * k
            cache["name"], cache["o"] = name, o
        return cache["o"]

    def truth_genes(name):
        """[(gid, start, end, strand, bases on read, gene length)] of the read's genome."""
        t = truth[name]
        out = []
        for gid, places in genes[t["acc"]].items():
            for s, e, strand in places:
                on = 0
                for rs, cs, ln in t["blocks"]:
                    if t["strand"] == "+":
                        lo, hi = cs, cs + ln
                    else:
                        lo, hi = cs - ln + 1, cs + 1
                    on += max(0, min(hi, e) - max(lo, s))
                if on > 0:
                    out.append((gid, s, e, strand, on, e - s))
        return out

    stats = collections.Counter()
    problems = collections.defaultdict(list)
    by_read = collections.defaultdict(list)
    mapq_ok, mapq_bad = collections.Counter(), collections.Counter()
    for r in parse_sam(a.sam):
        qname, flag, rname, pos, mapq, cigar, seq, qual = r[0], int(r[1]), r[2], int(r[3]), int(r[4]), r[5], r[9], r[10]
        tags = r[11]
        stats["records"] += 1
        if qname not in truth:
            stats["unknown_read"] += 1
            continue
        read, rqual = seqs[qname]
        L = len(read)
        ops = cigar_ops(cigar)
        clip_l = ops[0][0] if ops and ops[0][1] == "H" else 0
        clip_r = ops[-1][0] if len(ops) > 1 and ops[-1][1] == "H" else 0
        inner = ops[(1 if clip_l else 0):(len(ops) - 1 if clip_r else len(ops))]
        if any(op == "H" for _, op in inner):
            problems["H_inside"].append(r[:6])
        qlen = sum(n for n, op in inner if op in "MIS=X")
        if qlen != len(seq):
            problems["cigar_seq_len"].append(r[:6] + [len(seq)])
        if clip_l + len(seq) + clip_r != L:
            problems["clip_sum"].append(r[:6] + [clip_l, len(seq), clip_r, L])
        rev = bool(flag & 0x10)
        oriented = revcomp(read) if rev else read
        if seq != oriented[clip_l:clip_l + len(seq)]:
            problems["seq_mismatch"].append(r[:6])
        if a.fasta_qual:
            if set(qual) != {a.fasta_qual}:
                problems["fasta_qual"].append(r[:6] + [sorted(set(qual))[:5]])
        elif rqual is not None:
            oq = rqual[::-1] if rev else rqual
            if qual != oq[clip_l:clip_l + len(seq)]:
                problems["qual_mismatch"].append(r[:6])
        if flag & 0x4 or flag & 0x1:
            problems["bad_flag"].append(r[:6])
        tid, gid = (int(x) for x in rname.split("_"))
        t = truth[qname]
        true_tid = world.taxid[t["acc"]]
        secondary = bool(flag & 0x100)
        supp = bool(flag & 0x800)
        if secondary and supp:
            problems["secondary_and_supplementary"].append(r[:6])
        stats["secondary" if secondary else ("supplementary" if supp else "primary")] += 1
        if "ZR" in tags:
            stats["ZR"] += 1
            if secondary:
                problems["ZR_on_secondary"].append(r[:6])
        if secondary:
            if mapq != 0:
                problems["secondary_mapq"].append(r[:6])
            continue
        # coordinates
        origin = origin_of(qname)
        exact = near = total = 0
        places = genes[t["acc"]].get(gid, []) if tid == true_tid else []
        best_place = None
        if places:
            best = (-1, None)
            for place in places:
                gs, ge, gstrand = place
                ex = nr = tot = 0
                rp, gp = clip_l, pos - 1
                for n, op in inner:
                    if op in "M=X":
                        for k in range(n):
                            j = rp + k
                            i = (L - 1 - j) if rev else j
                            g = origin[i]
                            if g < 0:
                                continue
                            tot += 1
                            exp = gs + gp + k if gstrand == "+" else ge - 1 - (gp + k)
                            ex += g == exp
                            nr += abs(g - exp) <= 30
                        rp += n
                        gp += n
                    elif op in "IS":
                        rp += n
                    elif op == "D":
                        gp += n
                if nr > best[0]:
                    best = (nr, place, ex, tot)
            near, best_place, exact, total = best
        rep = world.rep[t["acc"]]
        correct = tid == true_tid and total > 0 and near >= 0.9 * total
        if "ZR" in tags:
            stats["ZR_correct" if correct else "ZR_wrong"] += 1
        if mapq < 4:
            stats["mapq_lt4_correct" if correct else "mapq_lt4_wrong"] += 1
        if correct:
            stats["records_correct_place"] += 1
            if rep:
                stats["rep_records"] += 1
                if exact < a.exact * total:
                    problems["rep_inexact_coordinates"].append(r[:6] + [exact, total])
            mapq_ok[min(mapq, 60) // 10 * 10] += 1
        else:
            stats["records_wrong"] += 1
            kind = "wrong_taxon" if tid != true_tid else ("gene_not_in_genome" if not places else "wrong_place")
            if kind == "wrong_taxon" and f"{true_tid}_{gid}" not in world.db_genes:
                kind = "wrong_taxon_own_lacks_gene"
            stats[kind] += 1
            problems[kind].append(r[:6] + [near, total, t["acc"]])
            mapq_bad[min(mapq, 60) // 10 * 10] += 1
        # reference span
        ref_span = sum(n for n, op in inner if op in "M=XD")
        by_read[qname].append((flag, tid, gid, best_place if correct else None, ref_span, mapq, "ZR" in tags))

    # per read
    for name, t in truth.items():
        recs = by_read.get(name, [])
        prim = [x for x in recs if not x[0] & 0x800]
        if recs and len(prim) != 1:
            problems["primaries_per_read"].append((name, len(prim)))
        stats["reads"] += 1
        stats["reads_with_records"] += bool(recs)
        true_tid = world.taxid[t["acc"]]
        for gid, s, e, strand, on, glen in truth_genes(name):
            full = on >= a.full * glen
            hits = [x for x in recs if x[1] == true_tid and x[2] == gid and x[3] is not None and x[3][0] == s]
            key = "full" if full else ("partial_ge50" if on >= 0.5 * glen else "partial_lt50")
            stats[key + "_genes"] += 1
            if hits:
                stats[key + "_found"] += 1
                if full and max(h[4] for h in hits) < 0.95 * glen:
                    problems["full_gene_short_alignment"].append((name, gid, glen, max(h[4] for h in hits)))
            elif full:
                if f"{true_tid}_{gid}" not in world.db_genes:
                    stats["full_missed_not_in_own_taxon_db"] += 1
                else:
                    problems["full_gene_missed"].append((name, gid, s, e, t["read_length"], t["strand"], t["start"]))
            if len(hits) > 1:
                problems["gene_twice"].append((name, gid, len(hits)))
    print("== stats")
    for k in sorted(stats):
        print(f"{k}\t{stats[k]}")
    if stats["full_genes"]:
        print(f"full_gene_recovery\t{stats['full_found'] / stats['full_genes']:.4f}")
    print("mapq (correct records, bins of 10, 60 = 60+)", dict(sorted(mapq_ok.items())))
    print("mapq (wrong records)", dict(sorted(mapq_bad.items())))
    print("== problems")
    for k, v in sorted(problems.items()):
        print(f"{k}\t{len(v)}")
        for x in v[:a.show]:
            print("   ", x)


if __name__ == "__main__":
    main()
