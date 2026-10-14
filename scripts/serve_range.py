#!/usr/bin/env python3
"""Static file server WITH byte-range support, for the offline video dashboards.

`python -m http.server` ignores Range headers, so a browser cannot seek an mp4 it serves
and the synchronised players stall. This adds single-range 206 responses and nothing else.

  serve_range.py <port> <directory>
"""
import http.server
import os
import re
import sys


class RangeHandler(http.server.SimpleHTTPRequestHandler):
    # Keep-alive: a page of many movies re-requests ranges constantly, and the browser's
    # per-host connection cap makes a fresh connection per range the bottleneck.
    protocol_version = 'HTTP/1.1'

    def send_head(self):
        rng = self.headers.get('Range')
        path = self.translate_path(self.path)
        if not rng or os.path.isdir(path) or not os.path.isfile(path):
            return super().send_head()
        m = re.fullmatch(r'bytes=(\d*)-(\d*)', rng.strip())
        size = os.path.getsize(path)
        if not m or (not m.group(1) and not m.group(2)):
            self.send_error(416)
            return None
        if m.group(1):
            start = int(m.group(1))
            end = min(int(m.group(2)), size-1) if m.group(2) else size-1
        else:                                   # suffix range: the last N bytes
            start, end = max(0, size-int(m.group(2))), size-1
        if start > end or start >= size:
            self.send_response(416)
            self.send_header('Content-Range', f'bytes */{size}')
            self.end_headers()
            return None
        f = open(path, 'rb')
        f.seek(start)
        self.send_response(206)
        self.send_header('Content-Type', self.guess_type(path))
        self.send_header('Content-Range', f'bytes {start}-{end}/{size}')
        self.send_header('Content-Length', str(end-start+1))
        self.send_header('Accept-Ranges', 'bytes')
        self.end_headers()
        self._remaining = end-start+1
        return f

    def copyfile(self, source, outputfile):
        left = getattr(self, '_remaining', None)
        if left is None:
            return super().copyfile(source, outputfile)
        self._remaining = None
        while left > 0:
            chunk = source.read(min(1 << 20, left))
            if not chunk:
                break
            try:
                outputfile.write(chunk)
            except (BrokenPipeError, ConnectionResetError):
                break                           # the browser cancelled the range; normal
            left -= len(chunk)

    def end_headers(self):
        if self.headers.get('Range') is None:
            self.send_header('Accept-Ranges', 'bytes')
        super().end_headers()


if __name__ == '__main__':
    port, root = int(sys.argv[1]), sys.argv[2]
    handler = lambda *a, **k: RangeHandler(*a, directory=root, **k)
    http.server.ThreadingHTTPServer(('127.0.0.1', port), handler).serve_forever()
