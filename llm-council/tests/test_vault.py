# ABOUTME: Tests for writing the council's synthesis as a note in the Obsidian vault.
# ABOUTME: Covers note naming, YAML frontmatter, collisions and a missing vault folder.

import json
import sys
from datetime import datetime
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.vault import (  # noqa: E402
    VaultError,
    build_frontmatter,
    note_path,
    note_title,
    write_vault_note,
)

NOW = datetime(2026, 10, 3, 15, 34, 12)


# --- naming ---------------------------------------------------------------


def test_note_title_keeps_a_plain_question_readable():
    assert note_title("Should we adopt a monorepo?") == "Should we adopt a monorepo"


def test_note_title_strips_characters_obsidian_cannot_link():
    title = note_title('a/b\\c:d*e?f"g<h>i|j#k^l[m]n')
    assert not set('/\\:*?"<>|#^[]') & set(title)
    assert title == "a b c d e f g h i j k l m n"


def test_note_title_uses_only_the_first_line_of_a_long_prompt():
    assert note_title("Summarise this\n\nlots of pasted text") == "Summarise this"


def test_note_title_is_bounded_and_breaks_on_a_word():
    title = note_title("word " * 60)
    assert len(title) <= 80
    assert not title.endswith(" ")
    assert title.split()[-1] == "word"


def test_note_title_falls_back_when_nothing_usable_remains():
    assert note_title("???") == "query"
    assert note_title(".hidden") == "hidden"


def test_note_path_leads_with_the_date_and_time(tmp_path):
    assert note_path("Why is the sky blue?", NOW, tmp_path) == (
        tmp_path / "2026-10-03 1534 Why is the sky blue.md"
    )


def test_note_path_never_overwrites_an_existing_note(tmp_path):
    (tmp_path / "2026-10-03 1534 q.md").write_text("earlier")
    (tmp_path / "2026-10-03 1534 q 2.md").write_text("earlier still")
    assert note_path("q", NOW, tmp_path) == tmp_path / "2026-10-03 1534 q 3.md"


# --- frontmatter ----------------------------------------------------------


def _frontmatter(**overrides):
    kwargs = dict(
        question='Is "X: Y" better?',
        now=NOW,
        panel=["claude claude-opus-5-5", "codex gpt-6-astra", "agy gemini-3.8-flash"],
        chair="claude claude-fable-5-1",
        ranking=["codex gpt-6-astra", "claude claude-opus-5-5", "agy gemini-3.8-flash"],
        mean_similarity=0.0268,
        run_dir=Path("/runs/20261003-153412-q"),
    )
    kwargs.update(overrides)
    return build_frontmatter(**kwargs)


def test_frontmatter_is_fenced_yaml_with_the_council_tag():
    fm = _frontmatter()
    assert fm.startswith("---\n") and fm.endswith("---\n")
    assert "tags:\n  - llm-council\n" in fm
    assert "created: 2026-10-03T15:34\n" in fm
    assert "type: LLM Council\n" in fm


def test_frontmatter_quotes_a_question_containing_yaml_syntax():
    fm = _frontmatter()
    assert f"question: {json.dumps('Is \"X: Y\" better?')}\n" in fm


def test_frontmatter_records_the_models_ranking_and_similarity():
    fm = _frontmatter()
    assert '  - "claude claude-opus-5-5"\n' in fm
    assert 'chair: "claude claude-fable-5-1"\n' in fm
    assert 'ranking:\n  - "codex gpt-6-astra"\n  - "claude claude-opus-5-5"\n' in fm
    assert "mean_similarity: 0.0268\n" in fm
    assert 'run_dir: "/runs/20261003-153412-q"\n' in fm


def test_frontmatter_records_an_empty_ranking_when_every_ballot_failed():
    assert "ranking: []\n" in _frontmatter(ranking=[])


# --- writing --------------------------------------------------------------


def _write(vault_dir, **overrides):
    kwargs = dict(
        final_text="# q\n\nThe answer.\n",
        question="q",
        run_dir=Path("/runs/r"),
        panel=["claude claude-opus-5-5"],
        chair="claude claude-fable-5-1",
        ranking=["claude claude-opus-5-5"],
        mean_similarity=0.05,
        vault_dir=vault_dir,
        now=NOW,
    )
    kwargs.update(overrides)
    return write_vault_note(**kwargs)


def test_write_vault_note_writes_frontmatter_then_the_final_document(tmp_path):
    path = _write(tmp_path)
    text = path.read_text()
    assert path == tmp_path / "2026-10-03 1534 q.md"
    assert text.startswith("---\n")
    body = text.split("---\n", 2)[2]
    assert body.startswith("\n# q\n\nThe answer.\n")


def test_write_vault_note_points_back_to_the_run_directory(tmp_path):
    text = _write(tmp_path).read_text()
    assert text.rstrip().endswith("- Run directory: `/runs/r`")


def test_write_vault_note_creates_the_council_folder_inside_an_existing_vault(tmp_path):
    path = _write(tmp_path / "LLM Council")
    assert path.parent == tmp_path / "LLM Council"
    assert path.exists()


def test_write_vault_note_refuses_to_invent_a_missing_vault(tmp_path):
    # An unmounted vault must fail loudly, not silently grow a stray directory tree.
    missing = tmp_path / "not-mounted" / "LLM Council"
    with pytest.raises(VaultError, match="not-mounted"):
        _write(missing)
    assert not (tmp_path / "not-mounted").exists()
