# ABOUTME: Tests for verify.sh, run against the real stand-in API over HTTP.
# ABOUTME: Covers the token checks and the check that a blocked machine gets no answer.

import os
import subprocess

from conftest import PROJECT, TEST_TOKEN, free_port

VERIFY = PROJECT / "verify.sh"


def run(*args, token=None):
    env = {key: value for key, value in os.environ.items() if key != "TOKEN"}
    if token is not None:
        env["TOKEN"] = token
    return subprocess.run(["bash", str(VERIFY), *args], env=env, capture_output=True, text=True, timeout=30)


def test_passes_when_the_api_enforces_the_token(api):
    result = run(api, token=TEST_TOKEN)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "without the token: 401" in result.stdout
    assert "with the token: 200" in result.stdout
    assert result.stderr == ""


def test_fails_when_the_token_is_wrong(api):
    result = run(api, token="not-the-token")
    assert result.returncode == 1
    assert "FAIL" in result.stdout
    assert "with the token: 401" in result.stdout


def test_never_prints_the_token(api):
    result = run(api, token=TEST_TOKEN)
    assert TEST_TOKEN not in result.stdout + result.stderr


def test_blocked_check_passes_when_nothing_answers():
    result = run("--blocked", f"http://127.0.0.1:{free_port()}")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "no answer" in result.stdout


def test_blocked_check_fails_when_the_api_answers(api):
    result = run("--blocked", api)
    assert result.returncode == 1
    assert "FAIL" in result.stdout


def test_token_check_needs_a_token(api):
    result = run(api)
    assert result.returncode == 2
    assert "TOKEN" in result.stderr


def test_needs_a_base_url():
    result = run(token=TEST_TOKEN)
    assert result.returncode == 2
    assert "Usage" in result.stderr
