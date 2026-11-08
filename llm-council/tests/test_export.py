# ABOUTME: Tests for exporting the chair's synthesis to markdown or docx.
# ABOUTME: Conversion shells out to pandoc via uv, so the project stays dependency-free.

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.export import (  # noqa: E402
    SUPPORTED_EXPORT_SUFFIXES,
    ExportError,
    export_command,
    export_document,
    validate_export_target,
)


class FakeConverter:
    def __init__(self, code=0, stderr=""):
        self.calls = []
        self.code = code
        self.stderr = stderr

    def __call__(self, cmd, timeout):
        self.calls.append(list(cmd))
        if self.code == 0:
            Path(cmd[-1]).parent.mkdir(parents=True, exist_ok=True)
            Path(cmd[-1]).write_bytes(b"PK fake docx")
        return self.code, "", self.stderr


# --- target validation ----------------------------------------------------


def test_markdown_and_docx_are_supported():
    assert {".md", ".docx"} <= SUPPORTED_EXPORT_SUFFIXES


def test_validate_rejects_an_unsupported_extension(tmp_path):
    with pytest.raises(ExportError) as excinfo:
        validate_export_target(tmp_path / "out.pdf")
    assert ".pdf" in str(excinfo.value)
    assert ".docx" in str(excinfo.value), "the error should name what is supported"


def test_validate_rejects_a_path_with_no_extension(tmp_path):
    with pytest.raises(ExportError):
        validate_export_target(tmp_path / "out")


def test_validate_accepts_a_target_whose_parent_does_not_exist_yet(tmp_path):
    """Parent directories are created at export time, like the run directory."""
    validate_export_target(tmp_path / "new" / "nested" / "out.docx")


def test_validate_rejects_an_existing_directory(tmp_path):
    with pytest.raises(ExportError):
        validate_export_target(tmp_path)


# --- markdown -------------------------------------------------------------


def test_markdown_export_copies_the_source_without_converting(tmp_path):
    source = tmp_path / "final.md"
    source.write_text("# Answer\n\nbody")
    target = tmp_path / "out" / "synthesis.md"
    converter = FakeConverter()

    export_document(source, target, runner=converter)

    assert target.read_text() == "# Answer\n\nbody"
    assert converter.calls == [], "markdown needs no converter"


def test_markdown_export_creates_missing_parent_directories(tmp_path):
    source = tmp_path / "final.md"
    source.write_text("x")
    target = tmp_path / "a" / "b" / "out.md"
    export_document(source, target, runner=FakeConverter())
    assert target.exists()


# --- docx -----------------------------------------------------------------


def test_export_command_runs_pandoc_through_uv_without_installing_it(tmp_path):
    """The project declares no dependency on pandoc; uv fetches it on demand."""
    cmd = export_command(tmp_path / "final.md", tmp_path / "out.docx")
    assert cmd[0] == "uv"
    assert "--with" in cmd
    assert "pypandoc_binary" in cmd
    assert str(tmp_path / "final.md") in cmd
    assert cmd[-1] == str(tmp_path / "out.docx")


def test_docx_export_invokes_the_converter(tmp_path):
    source = tmp_path / "final.md"
    source.write_text("# Answer")
    target = tmp_path / "synthesis.docx"
    converter = FakeConverter()

    export_document(source, target, runner=converter)

    assert len(converter.calls) == 1
    assert target.exists()


def test_docx_export_reports_a_conversion_failure(tmp_path):
    source = tmp_path / "final.md"
    source.write_text("# Answer")

    with pytest.raises(ExportError, match="no network"):
        export_document(
            source,
            tmp_path / "out.docx",
            runner=FakeConverter(code=1, stderr="no network"),
        )


def test_docx_export_errors_when_the_converter_produced_nothing(tmp_path):
    """A zero exit with no file would otherwise look like a successful export."""
    source = tmp_path / "final.md"
    source.write_text("# Answer")

    def silent_success(cmd, timeout):
        return 0, "", ""

    with pytest.raises(ExportError):
        export_document(source, tmp_path / "out.docx", runner=silent_success)


def test_missing_source_is_an_error(tmp_path):
    with pytest.raises(ExportError):
        export_document(tmp_path / "absent.md", tmp_path / "out.docx", runner=FakeConverter())
