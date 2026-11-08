# ABOUTME: Exports the chair's synthesis to a file outside the run directory.
# ABOUTME: Markdown is copied; docx is converted by pandoc fetched on demand through uv.

"""Export.

The project declares no dependency on pandoc. Conversion runs as
``uv run --with pypandoc_binary``, which fetches a wheel carrying the pandoc
binary into a temporary environment and caches it. That keeps the council
itself installable with nothing but a Python interpreter, at the cost of one
download the first time a docx is produced.

LibreOffice is not used: it has no markdown import filter and refuses to load a
``.md`` file at all.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Callable, Sequence

MARKDOWN_SUFFIXES = frozenset({".md", ".markdown"})
CONVERTED_SUFFIXES = frozenset({".docx"})
SUPPORTED_EXPORT_SUFFIXES = MARKDOWN_SUFFIXES | CONVERTED_SUFFIXES

Runner = Callable[[Sequence[str], int], "tuple[int, str, str]"]


class ExportError(RuntimeError):
    """The synthesis could not be written to the requested target."""


def _subprocess_runner(cmd: Sequence[str], timeout: int) -> tuple[int, str, str]:
    proc = subprocess.run(list(cmd), capture_output=True, text=True, timeout=timeout)
    return proc.returncode, proc.stdout, proc.stderr


def validate_export_target(target: Path | str) -> None:
    """Check the export path before the run starts.

    A council run takes minutes; discovering an unusable export path only at the
    end would waste all of it.
    """
    target = Path(target)
    supported = ", ".join(sorted(SUPPORTED_EXPORT_SUFFIXES))

    if target.is_dir():
        raise ExportError(f"export target is a directory: {target}")
    if target.suffix.lower() not in SUPPORTED_EXPORT_SUFFIXES:
        suffix = target.suffix or "(none)"
        raise ExportError(f"unsupported export format {suffix} for {target.name}; supported: {supported}")


def export_command(source: Path, target: Path) -> list[str]:
    """Build the pandoc-through-uv command.

    Paths are passed as argv to the child rather than interpolated into the
    inlined script, so a filename containing quotes cannot break it.
    """
    return [
        "uv",
        "run",
        "--with",
        "pypandoc_binary",
        "python",
        "-c",
        "import sys, pypandoc; pypandoc.convert_file(sys.argv[1], 'docx', outputfile=sys.argv[2])",
        str(source),
        str(target),
    ]


def export_document(
    source: Path | str,
    target: Path | str,
    runner: Runner = _subprocess_runner,
    timeout: int = 600,
) -> Path:
    source, target = Path(source), Path(target)
    validate_export_target(target)

    if not source.is_file():
        raise ExportError(f"nothing to export, {source} does not exist")

    target.parent.mkdir(parents=True, exist_ok=True)

    if target.suffix.lower() in MARKDOWN_SUFFIXES:
        shutil.copyfile(source, target)
        return target

    code, _, stderr = runner(export_command(source, target), timeout)
    if code != 0:
        raise ExportError(f"could not convert to {target.suffix}: {stderr.strip()}")
    if not target.exists():
        raise ExportError(f"conversion reported success but wrote nothing to {target}")
    return target
