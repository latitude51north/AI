# ABOUTME: Tests for firewall.sh, which limits a published port to one source in DOCKER-USER.
# ABOUTME: A recording stand-in replaces iptables, so no real firewall rule is ever changed.

import os
import subprocess

import pytest

from conftest import PROJECT

FIREWALL = PROJECT / "firewall.sh"

# Records every call. Remembers which of the two rules "exist", so -C, -I and -D
# behave as they would against a real chain.
FAKE_IPTABLES = """#!/usr/bin/env bash
echo "$*" >> "$FAKE_LOG"
kind=drop
[[ "$*" == *ACCEPT* ]] && kind=accept
case "$1" in
  -C) [[ -e "$FAKE_STATE/$kind" ]] ;;
  -I) touch "$FAKE_STATE/$kind" ;;
  -D) rm "$FAKE_STATE/$kind" ;;
  -L) echo "Chain DOCKER-USER" ;;
esac
"""

DROP = "DOCKER-USER -p tcp --dport 8081 -j DROP"
ACCEPT = "DOCKER-USER -p tcp --dport 8081 -s 192.0.2.20 -j ACCEPT"


@pytest.fixture
def chain(tmp_path):
    """Puts the recording iptables first on PATH. Returns helpers to run and inspect."""
    bin_dir = tmp_path / "bin"
    state = tmp_path / "state"
    bin_dir.mkdir()
    state.mkdir()
    fake = bin_dir / "iptables"
    fake.write_text(FAKE_IPTABLES)
    fake.chmod(0o755)
    log = tmp_path / "calls.log"
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}", "FAKE_LOG": str(log), "FAKE_STATE": str(state)}

    class Chain:
        def run(self, *args):
            return subprocess.run(["bash", str(FIREWALL), *args], env=env, capture_output=True, text=True, timeout=30)

        def calls(self):
            return log.read_text().splitlines() if log.exists() else []

        def forget_calls(self):
            log.unlink()

    return Chain()


def test_adds_drop_first_so_accept_ends_above_it(chain):
    result = chain.run("192.0.2.20")
    assert result.returncode == 0, result.stderr
    calls = chain.calls()
    assert calls.index(f"-I {DROP}") < calls.index(f"-I {ACCEPT}")
    assert calls[-1] == "-L DOCKER-USER -n --line-numbers"
    assert "Chain DOCKER-USER" in result.stdout


def test_running_twice_does_not_stack_duplicate_rules(chain):
    chain.run("192.0.2.20")
    chain.forget_calls()
    result = chain.run("192.0.2.20")
    assert result.returncode == 0, result.stderr
    calls = chain.calls()
    assert f"-I {DROP}" not in calls, "added a second DROP"
    assert calls.count(f"-D {ACCEPT}") == 1
    assert calls.count(f"-I {ACCEPT}") == 1
    assert calls.index(f"-D {ACCEPT}") < calls.index(f"-I {ACCEPT}")


def test_never_removes_the_drop_rule(chain):
    chain.run("192.0.2.20")
    chain.run("192.0.2.20")
    assert f"-I {DROP}" in chain.calls()
    assert not [call for call in chain.calls() if call.startswith("-D") and "DROP" in call]


def test_port_can_be_chosen(chain):
    result = chain.run("192.0.2.20", "9000")
    assert result.returncode == 0, result.stderr
    assert "-I DOCKER-USER -p tcp --dport 9000 -s 192.0.2.20 -j ACCEPT" in chain.calls()


def test_dry_run_prints_the_rules_and_changes_nothing(chain):
    result = chain.run("--dry-run", "192.0.2.20")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert lines.index(f"iptables -I {DROP}") < lines.index(f"iptables -I {ACCEPT}")
    assert chain.calls() == []


@pytest.mark.parametrize("source", ["not-an-address", "192.0.2", "192.0.2.300", "192.0.2.20/40", "192.0.2.20; reboot"])
def test_rejects_a_source_that_is_not_an_ipv4_address(chain, source):
    result = chain.run(source)
    assert result.returncode == 2
    assert "source" in result.stderr
    assert chain.calls() == []


@pytest.mark.parametrize("port", ["0", "70000", "http", "80a"])
def test_rejects_a_port_that_is_not_a_port(chain, port):
    result = chain.run("192.0.2.20", port)
    assert result.returncode == 2
    assert "port" in result.stderr
    assert chain.calls() == []


def test_accepts_a_network_as_the_source(chain):
    result = chain.run("192.0.2.0/24")
    assert result.returncode == 0, result.stderr
    assert "-I DOCKER-USER -p tcp --dport 8081 -s 192.0.2.0/24 -j ACCEPT" in chain.calls()


def test_needs_a_source(chain):
    result = chain.run()
    assert result.returncode == 2
    assert "Usage" in result.stderr
