"""Phase durations of a `protal --build` log timestamped by stamp.py, plus the STAT lines.
Usage: phases.py <label> <log>"""
import re
import sys

label, path = sys.argv[1:3]
marks = [("start", "Run Build"), ("pass1", "Build: iterate records"), ("pointers", "minimizers:"),
         ("pass2", "After first put"), ("unique", "Check Uniqueness"), ("stats", "Save unique kmer info"),
         ("write", "Write "), ("check", "Check"), ("end", "Run build took")]
seen, stats, uniq = {}, "", {"queries": 0, "scanned": 0, "single": 0}
for line in open(path):
    m = re.match(r"\s*([\d.]+) (.*)", line)
    if not m:
        continue
    t, text = float(m.group(1)), m.group(2)
    for name, prefix in marks:
        if name not in seen and (text == prefix if name == "check" else text.startswith(prefix)):
            seen[name] = t
    if text.startswith("STAT unique_kmers"):
        stats = text[len("STAT unique_kmers "):]
    if text.startswith("STAT uniqueness thread"):
        for k in uniq:
            uniq[k] += int(re.search(rf"{k}=(\d+)", text).group(1))
names = [n for n, _ in marks]
out = [label, f"init={seen.get('start', 0):.2f}"]
for a, b in zip(names, names[1:]):
    out.append(f"{a}={seen[b] - seen[a]:.2f}" if a in seen and b in seen else f"{a}=?")
out.append("total=" + (f"{seen['end']:.2f}" if "end" in seen else "?"))
print(" ".join(out))
print("   ", stats, " | uniqueness", " ".join(f"{k}={v}" for k, v in uniq.items()))
