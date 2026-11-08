# ABOUTME: Shared pytest fixtures for the council test suite.
# ABOUTME: Redirects the Obsidian vault so no test can write into the real one.

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import council.vault  # noqa: E402


@pytest.fixture(autouse=True)
def isolated_vault(tmp_path, monkeypatch):
    """Point every run at a throwaway vault folder unless a test passes its own."""
    vault_dir = tmp_path / "vault"
    monkeypatch.setattr(council.vault, "VAULT_DIR", vault_dir)
    return vault_dir
