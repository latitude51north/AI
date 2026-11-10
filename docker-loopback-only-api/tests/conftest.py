# ABOUTME: Shared fixtures for the loopback-only API example's tests.
# ABOUTME: Starts the real stand-in API on a free loopback port for tests that need a live server.

import os
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest

PROJECT = Path(__file__).resolve().parent.parent

# A made-up value used only by the tests. It protects nothing.
TEST_TOKEN = "token-used-only-by-the-tests"


def free_port():
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        return probe.getsockname()[1]


def wait_until_listening(port, process, timeout=10):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError(f"the server exited early with code {process.returncode}")
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.2):
                return
        except OSError:
            time.sleep(0.05)
    raise RuntimeError("the server did not start listening in time")


@pytest.fixture
def project():
    return PROJECT


@pytest.fixture
def api():
    """The stand-in API running as its own process. Yields its base URL."""
    port = free_port()
    process = subprocess.Popen(
        [sys.executable, str(PROJECT / "app" / "server.py")],
        env={**os.environ, "API_TOKEN": TEST_TOKEN, "PORT": str(port)},
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        wait_until_listening(port, process)
        yield f"http://127.0.0.1:{port}"
    finally:
        process.terminate()
        process.wait(timeout=10)
