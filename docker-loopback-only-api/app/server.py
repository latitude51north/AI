#!/usr/bin/env python3
# ABOUTME: Stand-in for a service that listens only on loopback and requires a bearer token.
# ABOUTME: It plays the app's API in this example. Replace it with the real service.

"""A small HTTP API that can only be reached on 127.0.0.1.

It behaves like the kind of service the example is about: it binds to loopback,
there is no setting to change that, and every request must carry a bearer token.

Settings come from the environment:

    API_TOKEN   the bearer token callers must send. Required.
    PORT        the loopback port to listen on. Defaults to 8080.
"""

import hmac
import json
import os
import signal
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

# Loopback only, on purpose. This is the constraint the sidecar works around.
HOST = "127.0.0.1"
DEFAULT_PORT = 8080


def make_server(token, port=DEFAULT_PORT):
    """Return a server on loopback that answers only callers who send the token."""
    expected = f"Bearer {token}".encode()

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            supplied = self.headers.get("Authorization", "").encode()
            if not hmac.compare_digest(supplied, expected):
                self.reply(401, {"error": "unauthorised"}, {"WWW-Authenticate": "Bearer"})
            elif self.path == "/api/health":
                self.reply(200, {"status": "ok"})
            else:
                self.reply(404, {"error": "not found"})

        def reply(self, status, payload, extra_headers=None):
            body = json.dumps(payload).encode()
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            for name, value in (extra_headers or {}).items():
                self.send_header(name, value)
            self.end_headers()
            self.wfile.write(body)

    return ThreadingHTTPServer((HOST, port), Handler)


def stop(signum, frame):
    sys.exit(0)


def main():
    token = os.environ.get("API_TOKEN")
    if not token:
        sys.exit("API_TOKEN is not set. The API will not start without a token.")
    port = int(os.environ.get("PORT", DEFAULT_PORT))

    # As process 1 in a container, Python ignores SIGTERM unless it is handled,
    # and "docker stop" would wait out its grace period before killing the process.
    signal.signal(signal.SIGTERM, stop)

    with make_server(token, port) as server:
        print(f"Listening on http://{HOST}:{port}", file=sys.stderr, flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
