# ABOUTME: Loads local files the council should reason over, converting them to markdown.
# ABOUTME: Markdown is read directly; documents go through `firecrawl parse`.

"""Attachments.

An attachment is treated as another cited source, tagged ``[A1]``, ``[A2]`` and so
on alongside the web brief's ``[S1]``. The converted text is inlined into every
member's prompt rather than left on disk for members to open, which keeps access
symmetric: no member needs file tools, so none has an advantage over another.

Attachments are also shown to reviewers in stage 2. Without that, a reviewer
could not check whether a claim about an attached document is accurate — the
same reason the evidence brief is repeated at review time.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

from council.research import TRUNCATION_MARKER, slugify, truncate_markdown

# Read straight off disk: Firecrawl cannot parse these, and spending credits to
# convert text that is already text would be waste.
TEXT_SUFFIXES = frozenset({".md", ".markdown", ".txt"})

# Converted by `firecrawl parse`, which covers these document formats.
PARSE_SUFFIXES = frozenset(
    {".pdf", ".docx", ".doc", ".odt", ".rtf", ".xlsx", ".xls", ".html", ".htm"}
)

SUPPORTED_SUFFIXES = TEXT_SUFFIXES | PARSE_SUFFIXES

# Generous per-document budget. Gemini's whole prompt is capped at 96 KB by
# MAX_ARG_STRLEN, so several large attachments plus a brief will still raise
# PromptTooLargeError at dispatch; this only stops one document eating it all.
DEFAULT_MAX_CHARS_PER_ATTACHMENT = 20_000

Runner = Callable[[Sequence[str], int], "tuple[int, str, str]"]


class AttachmentError(RuntimeError):
    """An attachment could not be read, converted, or is not a supported type."""


@dataclass(frozen=True)
class Attachment:
    label: str
    title: str
    source_path: Path
    markdown: str


def _subprocess_runner(cmd: Sequence[str], timeout: int) -> tuple[int, str, str]:
    proc = subprocess.run(list(cmd), capture_output=True, text=True, timeout=timeout)
    return proc.returncode, proc.stdout, proc.stderr


def validate_attachment_paths(paths: Sequence[Path | str]) -> None:
    """Check every attachment up front.

    Runs before retrieval so a typo in a filename fails immediately rather than
    after a run has already spent Firecrawl credits and subscription quota.
    """
    supported = ", ".join(sorted(SUPPORTED_SUFFIXES))
    for raw in paths:
        path = Path(raw)
        if not path.exists():
            raise AttachmentError(f"attachment does not exist: {path}")
        if not path.is_file():
            raise AttachmentError(f"attachment is not a file: {path}")
        if path.suffix.lower() not in SUPPORTED_SUFFIXES:
            raise AttachmentError(
                f"unsupported attachment type {path.suffix!r} for {path.name}; "
                f"supported: {supported}"
            )


def parse_command(source: Path, out_path: Path) -> list[str]:
    return [
        "firecrawl",
        "parse",
        str(source),
        "--format",
        "markdown",
        "--only-main-content",
        "-o",
        str(out_path),
    ]


def load_attachments(
    paths: Sequence[Path | str],
    out_dir: Path | str,
    runner: Runner = _subprocess_runner,
    timeout: int = 300,
) -> list[Attachment]:
    """Convert each attachment to markdown and record it in the run directory."""
    paths = [Path(p) for p in paths]
    validate_attachment_paths(paths)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    attachments: list[Attachment] = []
    for index, source in enumerate(paths, start=1):
        label = f"A{index}"
        target = out_dir / f"{label}-{slugify(source.stem, maxlen=48)}.md"

        if source.suffix.lower() in TEXT_SUFFIXES:
            text = source.read_text(errors="replace")
        else:
            code, _, stderr = runner(parse_command(source, target), timeout)
            if code != 0:
                raise AttachmentError(f"could not parse {source.name}: {stderr.strip()}")
            text = target.read_text(errors="replace") if target.exists() else ""

        if not text.strip():
            raise AttachmentError(
                f"{source.name} produced no text; it may be scanned images or empty"
            )

        target.write_text(text)
        attachments.append(
            Attachment(
                label=label,
                title=source.name,
                source_path=source,
                markdown=text.strip(),
            )
        )
    return attachments


def build_attachment_block(
    attachments: Sequence[Attachment],
    max_chars_each: int = DEFAULT_MAX_CHARS_PER_ATTACHMENT,
) -> str:
    if not attachments:
        return ""
    blocks = []
    for attachment in attachments:
        body, _ = truncate_markdown(attachment.markdown, max_chars_each)
        blocks.append(f"[{attachment.label}] {attachment.title}\n\n{body}")
    return "\n\n---\n\n".join(blocks)


__all__ = [
    "Attachment",
    "AttachmentError",
    "DEFAULT_MAX_CHARS_PER_ATTACHMENT",
    "PARSE_SUFFIXES",
    "SUPPORTED_SUFFIXES",
    "TEXT_SUFFIXES",
    "TRUNCATION_MARKER",
    "build_attachment_block",
    "load_attachments",
    "parse_command",
    "validate_attachment_paths",
]
