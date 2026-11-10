# ABOUTME: Tests for the stand-in API: where it listens and how it treats the bearer token.
# ABOUTME: Requests go over real HTTP to a real server process.

import importlib.util
import json
import os
import subprocess
import sys
import urllib.error
import urllib.request

from conftest import PROJECT, TEST_TOKEN, free_port, wait_until_listening

SERVER = PROJECT / "app" / "server.py"


def load_server_module():
    spec = importlib.util.spec_from_file_location("standin_server", SERVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def get(url, token=None):
    request = urllib.request.Request(url)
    if token is not None:
        request.add_header("Authorization", f"Bearer {token}")
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, dict(response.headers), response.read()
    except urllib.error.HTTPError as error:
        return error.code, dict(error.headers), error.read()


def test_binds_loopback_only_on_the_port_the_proxy_expects():
    module = load_server_module()
    assert module.DEFAULT_PORT == 8080
    with module.make_server(TEST_TOKEN, port=0) as server:
        assert server.server_address[0] == "127.0.0.1"


def test_request_without_a_token_is_refused(api):
    status, headers, _ = get(f"{api}/api/health")
    assert status == 401
    assert headers["WWW-Authenticate"] == "Bearer"


def test_request_with_the_wrong_token_is_refused(api):
    status, _, _ = get(f"{api}/api/health", token="not-the-token")
    assert status == 401


def test_request_with_the_token_is_answered(api):
    status, headers, body = get(f"{api}/api/health", token=TEST_TOKEN)
    assert status == 200
    assert headers["Content-Type"] == "application/json"
    assert json.loads(body) == {"status": "ok"}


def test_unknown_path_is_not_found(api):
    status, _, _ = get(f"{api}/api/nothing-here", token=TEST_TOKEN)
    assert status == 404


def test_unknown_path_without_a_token_is_still_refused(api):
    status, _, _ = get(f"{api}/api/nothing-here")
    assert status == 401


def test_refuses_to_start_without_a_token():
    env = {key: value for key, value in os.environ.items() if key != "API_TOKEN"}
    result = subprocess.run([sys.executable, str(SERVER)], env=env, capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert "API_TOKEN" in result.stderr


def test_stops_cleanly_when_asked_to_terminate():
    # In a container the server is process 1, which ignores SIGTERM unless it handles it.
    port = free_port()
    process = subprocess.Popen(
        [sys.executable, str(SERVER)],
        env={**os.environ, "API_TOKEN": TEST_TOKEN, "PORT": str(port)},
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        wait_until_listening(port, process)
        process.terminate()
        assert process.wait(timeout=10) == 0
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
