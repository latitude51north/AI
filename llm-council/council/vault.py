# ABOUTME: Writes the chair's synthesis as a note in the personal Obsidian vault.
# ABOUTME: Adds YAML frontmatter for search and Dataview, and never overwrites a note.

"""Vault notes.

Every run leaves one note in the vault so answers are read where the rest of
the notes live. Only the synthesis goes there; the stage files, prompts and
JSON stay in the run directory, which the note links back to.

The vault folder is created if absent, but the vault itself is not: if the
vault is not mounted the write fails rather than quietly building a stray
directory tree at the vault's path.
"""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path
from typing import Sequence

VAULT_DIR = Path("/config/Documents/Harvinder remote vault/copilot/LLM Council")

# Characters Obsidian forbids in a file name or that break [[wikilinks]].
UNSAFE_TITLE_CHARS = re.compile(r'[\\/:*?"<>|#^\[\]]')
MAX_TITLE_LENGTH = 80


class VaultError(RuntimeError):
    """The note could not be written to the vault."""


def note_title(question: str) -> str:
    """A readable file-name title from the question's first line."""
    first_line = question.strip().splitlines()[0] if question.strip() else ""
    title = " ".join(UNSAFE_TITLE_CHARS.sub(" ", first_line).split()).lstrip(".").strip()
    if len(title) > MAX_TITLE_LENGTH:
        title = title[: MAX_TITLE_LENGTH + 1].rsplit(" ", 1)[0].rstrip()
    return title or "query"


def note_path(question: str, now: datetime, vault_dir: Path) -> Path:
    """The next free path for this question, suffixed rather than overwriting."""
    stem = f"{now:%Y-%m-%d %H%M} {note_title(question)}"
    path = vault_dir / f"{stem}.md"
    n = 2
    while path.exists():
        path = vault_dir / f"{stem} {n}.md"
        n += 1
    return path


def _yaml_list(key: str, values: Sequence[str]) -> list[str]:
    if not values:
        return [f"{key}: []"]
    return [f"{key}:"] + [f"  - {json.dumps(v)}" for v in values]


def build_frontmatter(
    question: str,
    now: datetime,
    panel: Sequence[str],
    chair: str,
    ranking: Sequence[str],
    mean_similarity: float | None,
    run_dir: Path,
) -> str:
    """YAML frontmatter. Strings are JSON-quoted, which YAML reads as double-quoted scalars."""
    lines = [
        "---",
        f"created: {now:%Y-%m-%dT%H:%M}",
        "tags:",
        "  - llm-council",
        "type: LLM Council",
        f"question: {json.dumps(question)}",
        *_yaml_list("panel", panel),
        f"chair: {json.dumps(chair)}",
        *_yaml_list("ranking", ranking),
        f"mean_similarity: {json.dumps(mean_similarity)}",
        f"run_dir: {json.dumps(str(run_dir))}",
        "---",
    ]
    return "\n".join(lines) + "\n"


def write_vault_note(
    final_text: str,
    question: str,
    run_dir: Path,
    panel: Sequence[str],
    chair: str,
    ranking: Sequence[str],
    mean_similarity: float | None,
    vault_dir: Path = VAULT_DIR,
    now: datetime | None = None,
) -> Path:
    now = now or datetime.now()
    vault_dir = Path(vault_dir)
    try:
        vault_dir.mkdir(exist_ok=True)
    except FileNotFoundError as exc:
        raise VaultError(f"vault is not available at {vault_dir.parent}") from exc

    path = note_path(question, now, vault_dir)
    frontmatter = build_frontmatter(
        question, now, panel, chair, ranking, mean_similarity, run_dir
    )
    body = final_text.rstrip() + f"\n- Run directory: `{run_dir}`\n"
    path.write_text(frontmatter + "\n" + body)
    return path
