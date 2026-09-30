#!/usr/bin/env python3
"""A localhost HTTP server for download_gtdb.py tests: serves a folder with Range support (206, and 416
for a range starting at or after the end, as nginx/Apache do), and fault injection read from FAULTS
(a JSON file re-read per request): {"truncate": {"<path suffix>": bytes}} sends a Content-Length of the
whole file but closes after that many bytes; {"no_range": true} ignores Range (200)."""
import http.server
import json
import os
import sys
import threading

ROOT, PORT, FAULTS = sys.argv[1], int(sys.argv[2]), sys.argv[3]


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=ROOT, **k)

    def log_message(self, fmt, *args):
        sys.stderr.write("SERVER " + (fmt % args) + " range=" + str(self.headers.get("Range")) + "\n")

    def do_GET(self):
        path = self.translate_path(self.path)
        if os.path.isdir(path):
            return super().do_GET()
        try:
            faults = json.load(open(FAULTS))
        except Exception:
            faults = {}
        if not os.path.isfile(path):
            self.send_error(404)
            return
        size = os.path.getsize(path)
        start = 0
        rng = self.headers.get("Range")
        if rng and not faults.get("no_range"):
            start = int(rng.split("=")[1].split("-")[0])
            if start >= size:
                self.send_response(416)
                self.send_header("Content-Range", f"bytes */{size}")
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{size - 1}/{size}")
        else:
            self.send_response(200)
        self.send_header("Content-Length", str(size - start))
        self.end_headers()
        limit = None
        for suffix, n in faults.get("truncate", {}).items():
            if self.path.endswith(suffix):
                limit = n
        with open(path, "rb") as fh:
            fh.seek(start)
            data = fh.read()
        if limit is not None:
            data = data[:limit]
        self.wfile.write(data)
        self.wfile.flush()
        if limit is not None:
            self.close_connection = True


server = http.server.ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
server.serve_forever()
