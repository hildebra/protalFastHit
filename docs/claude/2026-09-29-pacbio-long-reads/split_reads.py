#!/usr/bin/env python3
"""split_reads.py - cut every read of a FASTQ into consecutive pieces of --length bases.

A piece is named <read>_p<i> (the read it came from stays recoverable); a last piece shorter than
--min is dropped. This is the "split per default" alternative to aligning long reads whole.
"""

import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--min", type=int, default=50)
    args = ap.parse_args()
    pieces = 0
    with open(args.fastq) as fh, open(args.out, "w") as out:
        while True:
            header = fh.readline().rstrip("\n")
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            name = header[1:].split()[0]
            for i, start in enumerate(range(0, len(seq), args.length)):
                piece = seq[start:start + args.length]
                if len(piece) < args.min:
                    continue
                out.write(f"@{name}_p{i}\n{piece}\n+\n{qual[start:start + args.length]}\n")
                pieces += 1
    print(f"{pieces} pieces of {args.length} bp")


if __name__ == "__main__":
    main()
