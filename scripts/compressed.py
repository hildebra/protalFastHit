#!/usr/bin/env python3
"""Compressed files for the training data collection's reads and templates: zstd or gzip.

open_write(path) writes zstd if the name ends in .zst, gzip if in .gz (a trailing .partial does not count), plain
otherwise; open_read(path) reads a file by its first bytes (zstd, gzip or plain), as protal reads read files. zstd is
Python 3.14's compression.zstd where there is one, else the zstd command through a pipe (in the protal-db-build
environment and on Ubuntu: apt install zstd).

usage: import compressed; with compressed.open_write("reads.fq.zst") as fh: fh.write(b"...")
"""

import gzip
import shutil
import subprocess

ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"
SUFFIXES = {"zstd": ".zst", "gzip": ".gz"}
LEVELS = {"zstd": 3, "gzip": 1}  # fast: the files are read once, by protal


def compression_of(path):
    """zstd, gzip or None by a file name (without a trailing .partial)."""
    name = path[:-len(".partial")] if path.endswith(".partial") else path
    return "zstd" if name.endswith(".zst") else "gzip" if name.endswith(".gz") else None


def is_zstd(head):
    """Whether bytes start a zstd frame or a skippable frame."""
    return head[:4] == ZSTD_MAGIC or (len(head) >= 4 and head[0] & 0xf0 == 0x50 and head[1:4] == b"\x2a\x4d\x18")


def _stdlib_zstd():
    try:
        from compression import zstd  # Python 3.14
        return zstd
    except ImportError:
        return None


def _zstd_command():
    command = shutil.which("zstd")
    if not command:
        raise OSError("zstd files need Python 3.14 (compression.zstd) or the zstd command on PATH")
    return command


class _Pipe:
    """A binary file object over a zstd process: written to (compressing into a file) or read from (decompressing
    one). close() waits for the process and raises OSError if it failed."""

    def __init__(self, command, writing):
        self.writing = writing
        self.process = subprocess.Popen(command, stdin=subprocess.PIPE if writing else subprocess.DEVNULL,
                                        stdout=None if writing else subprocess.PIPE, stderr=subprocess.PIPE)
        self.stream = self.process.stdin if writing else self.process.stdout
        self.command = command

    def write(self, data):
        return self.stream.write(data)

    def read(self, size=-1):
        return self.stream.read(size)

    def readline(self, size=-1):
        return self.stream.readline(size)

    def __iter__(self):
        return iter(self.stream)

    def close(self):
        if self.process is None:
            return
        process, self.process = self.process, None
        if self.writing:
            process.stdin.close()
        else:
            process.stdout.read()  # to the end, so that the process does not wait on a full pipe
            process.stdout.close()
        error = process.stderr.read().decode(errors="replace").strip()
        process.stderr.close()
        if process.wait() != 0:
            raise OSError(f"{' '.join(self.command)} failed ({process.returncode}): {error[-300:]}")

    def __enter__(self):
        return self

    def __exit__(self, kind, value, traceback):
        if kind is None:
            self.close()
        else:  # an exception already: the process ends without hiding it
            try:
                self.close()
            except OSError:
                pass


def open_write(path, level=None):
    """A binary file object writing `path`: zstd (.zst), gzip (.gz) or plain, by its name; level: the codec's
    (default LEVELS)."""
    kind = compression_of(path)
    if kind is None:
        return open(path, "wb")
    level = LEVELS[kind] if level is None else level
    if kind == "gzip":
        return gzip.open(path, "wb", compresslevel=level)
    zstd = _stdlib_zstd()
    if zstd is not None:
        return zstd.open(path, "wb", level=level)
    return _Pipe([_zstd_command(), "-q", "-f", f"-{level}", "-o", path], writing=True)


def open_read(path):
    """A binary file object reading `path` decompressed: zstd, gzip or plain, by its first bytes."""
    with open(path, "rb") as fh:
        head = fh.read(4)
    if head[:2] == b"\x1f\x8b":
        return gzip.open(path, "rb")
    if is_zstd(head):
        zstd = _stdlib_zstd()
        if zstd is not None:
            return zstd.open(path, "rb")
        return _Pipe([_zstd_command(), "-dcq", path], writing=False)
    return open(path, "rb")


def read_text(path):
    """The decompressed content of a file, as text."""
    with open_read(path) as fh:
        return fh.read().decode()
