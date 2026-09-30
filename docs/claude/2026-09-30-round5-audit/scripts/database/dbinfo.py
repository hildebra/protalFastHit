#!/usr/bin/env python3
"""List the members of a database.protal: seek table, directory frame, each member's frames.

Read-only inspection for the audit; decompresses frame 0 with the zstd CLI.
Usage: dbinfo.py FILE [--frames]
"""
import argparse
import struct
import subprocess


def seek_table(path):
    with open(path, "rb") as f:
        data = f.read()
    n, desc, magic = struct.unpack("<IBI", data[-9:])
    assert magic == 0x8F92EAB1, "no seek table"
    entry = 12 if desc & 0x80 else 8
    size = n * entry + 9
    start = len(data) - size - 8
    frames, off = [], 0
    for i in range(n):
        c, d = struct.unpack("<II", data[start + 8 + i * entry:start + 16 + i * entry])
        frames.append((off, c, d))
        off += c
    return data, frames


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file")
    ap.add_argument("--frames", action="store_true")
    a = ap.parse_args()
    data, frames = seek_table(a.file)
    off, c, d = frames[0]
    directory = subprocess.run(["zstd", "-dc"], input=data[off:off + c], capture_output=True, check=True).stdout
    assert directory[:8] == b"PROTALDB"
    version, count = struct.unpack("<QQ", directory[8:24])
    pos = 24
    print(f"frames={len(frames)} version={version} members={count}")
    for _ in range(count):
        (length,) = struct.unpack("<Q", directory[pos:pos + 8]); pos += 8
        name = directory[pos:pos + length].decode(); pos += length
        first, nf = struct.unpack("<QQ", directory[pos:pos + 16]); pos += 16
        sub = frames[first:first + nf]
        start = sub[0][0] if sub else 0
        print(f"{name}\tfirst={first}\tframes={nf}\tfile_offset={start}\tcompressed={sum(x[1] for x in sub)}\tsize={sum(x[2] for x in sub)}")
        if a.frames:
            for i, (o, cc, dd) in enumerate(sub):
                print(f"   frame {first + i}: offset={o} compressed={cc} content={dd}")


if __name__ == "__main__":
    main()
