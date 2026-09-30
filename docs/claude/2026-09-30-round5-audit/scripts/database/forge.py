#!/usr/bin/env python3
"""Forge a database.protal: keep the original data frames, replace the directory (frame 0) and/or
the seek table, to test Bundle::Open's validation. Read-only on the source; writes a new file.

Subcommands operate on the mini bundle's 6 members (index.prx, reference.fna, reference.map,
internal_taxonomy.dmp, unique_kmers.tsv, model_pe.xml)."""
import struct, subprocess, sys

SEEK_MAGIC = 0x8F92EAB1
FRAME_MAGIC = 0x184D2A5E


def parse(path):
    data = open(path, "rb").read()
    n, desc, magic = struct.unpack("<IBI", data[-9:])
    assert magic == SEEK_MAGIC
    entry = 12 if desc & 0x80 else 8
    size = n * entry + 9
    start = len(data) - size - 8
    frames, off = [], 0
    for i in range(n):
        c, d = struct.unpack("<II", data[start + 8 + i * entry:start + 16 + i * entry])
        frames.append([off, c, d])   # file_offset, compressed, content
        off += c
    return data, frames


def members(data, frames):
    o, c, d = frames[0]
    directory = subprocess.run(["zstd", "-dc"], input=data[o:o + c], capture_output=True, check=True).stdout
    assert directory[:8] == b"PROTALDB"
    version, count = struct.unpack("<QQ", directory[8:24]); pos = 24
    mem = []
    for _ in range(count):
        (ln,) = struct.unpack("<Q", directory[pos:pos + 8]); pos += 8
        name = directory[pos:pos + ln]; pos += ln
        first, nf = struct.unpack("<QQ", directory[pos:pos + 16]); pos += 16
        mem.append([name, first, nf])
    return version, mem


def build_directory(version, mem):
    out = bytearray(b"PROTALDB")
    out += struct.pack("<QQ", version, len(mem))
    for name, first, nf in mem:
        out += struct.pack("<Q", len(name)) + name + struct.pack("<QQ", first, nf)
    return bytes(out)


def emit(path, data, frames, new_dir=None, seek_overrides=None):
    """Re-emit with an optional replaced directory and optional seek-table content-size overrides."""
    frame_bytes = []
    for o, c, d in frames:
        frame_bytes.append(bytearray(data[o:o + c]))
    contents = [d for _, _, d in frames]
    if new_dir is not None:
        comp = subprocess.run(["zstd", "-c", "-3"], input=new_dir, capture_output=True, check=True).stdout
        frame_bytes[0] = bytearray(comp)
        contents[0] = len(new_dir)
    if seek_overrides:
        for idx, val in seek_overrides.items():
            contents[idx] = val
    body = b"".join(bytes(b) for b in frame_bytes)
    seek = bytearray()
    seek += struct.pack("<I", FRAME_MAGIC)
    table = bytearray()
    for b, d in zip(frame_bytes, contents):
        table += struct.pack("<II", len(b), d)
    seek += struct.pack("<I", len(table) + 9)
    seek += table
    seek += struct.pack("<I", len(frame_bytes))
    seek += b"\x00"
    seek += struct.pack("<I", SEEK_MAGIC)
    open(path, "wb").write(body + bytes(seek))


def main():
    src, cmd, dst = sys.argv[1], sys.argv[2], sys.argv[3]
    data, frames = parse(src)
    version, mem = members(data, frames)
    if cmd == "traversal":
        mem[2][0] = b"../escapee.map"           # reference.map -> path traversal name
        emit(dst, data, frames, build_directory(version, mem))
    elif cmd == "absolute":
        mem[2][0] = b"/tmp/escapee.map"
        emit(dst, data, frames, build_directory(version, mem))
    elif cmd == "dupe":
        mem[3][0] = mem[2][0]                    # two members with the same name
        emit(dst, data, frames, build_directory(version, mem))
    elif cmd == "notile":
        mem[2][1] += 1                           # first != previous next
        emit(dst, data, frames, build_directory(version, mem))
    elif cmd == "badcount":
        version_out = version
        mem2 = mem + [[b"extra.xml", frames.__len__(), 1]]  # claims a member past the frames
        emit(dst, data, frames, build_directory(version_out, mem2))
    elif cmd == "version":
        emit(dst, data, frames, build_directory(999, mem))
    elif cmd == "hugecontent":
        # Forge a data frame's decompressed size in the seek table to ~4 GB (index chunk frame 5).
        emit(dst, data, frames, seek_overrides={5: 0xFFFFFFF0})
    elif cmd == "hugedir":
        # Forge frame 0's (directory) content size to > 16 MB; Bundle::Open should reject pre-alloc.
        emit(dst, data, frames, seek_overrides={0: 0x11000000})
    else:
        print("unknown", cmd); sys.exit(2)
    print("wrote", dst)


if __name__ == "__main__":
    main()
