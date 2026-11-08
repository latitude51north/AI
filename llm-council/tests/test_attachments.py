# ABOUTME: Tests for local file attachments passed to the council.
# ABOUTME: Covers suffix routing, fail-fast validation, parsing and prompt assembly.

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.attachments import (  # noqa: E402
    PARSE_SUFFIXES,
    TEXT_SUFFIXES,
    Attachment,
    AttachmentError,
    build_attachment_block,
    load_attachments,
    parse_command,
    validate_attachment_paths,
)


class FakeParser:
    """Stands in for `firecrawl parse`, recording what it was asked to convert."""

    def __init__(self, output="# parsed\n\nbody from the document"):
        self.calls = []
        self.output = output

    def __call__(self, cmd, timeout):
        self.calls.append(list(cmd))
        out_path = Path(cmd[cmd.index("-o") + 1])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(self.output)
        return 0, "", ""


# --- validation (fail fast) -----------------------------------------------


def test_validate_rejects_a_missing_file(tmp_path):
    """Fail before the run starts rather than after spending credits on retrieval."""
    with pytest.raises(AttachmentError, match="does not exist"):
        validate_attachment_paths([tmp_path / "nope.pdf"])


def test_validate_rejects_an_unsupported_type(tmp_path):
    bad = tmp_path / "clip.png"
    bad.write_text("x")
    with pytest.raises(AttachmentError) as excinfo:
        validate_attachment_paths([bad])
    assert ".png" in str(excinfo.value)


def test_validate_names_every_supported_suffix_in_its_error(tmp_path):
    bad = tmp_path / "clip.png"
    bad.write_text("x")
    with pytest.raises(AttachmentError) as excinfo:
        validate_attachment_paths([bad])
    assert ".pdf" in str(excinfo.value) and ".md" in str(excinfo.value)


def test_validate_accepts_the_required_formats(tmp_path):
    paths = []
    for name in ("notes.md", "report.pdf", "spec.docx"):
        p = tmp_path / name
        p.write_text("x")
        paths.append(p)
    validate_attachment_paths(paths)


def test_validate_rejects_a_directory(tmp_path):
    with pytest.raises(AttachmentError):
        validate_attachment_paths([tmp_path])


# --- suffix routing -------------------------------------------------------


def test_markdown_is_read_directly_without_calling_firecrawl(tmp_path):
    """Firecrawl cannot parse .md and charging credits for a local read would be waste."""
    src = tmp_path / "notes.md"
    src.write_text("# Local notes\n\nalready markdown")
    parser = FakeParser()

    attachments = load_attachments([src], out_dir=tmp_path / "out", runner=parser)

    assert parser.calls == []
    assert "already markdown" in attachments[0].markdown


def test_pdf_is_routed_through_firecrawl_parse(tmp_path):
    src = tmp_path / "report.pdf"
    src.write_bytes(b"%PDF-1.4 fake")
    parser = FakeParser()

    attachments = load_attachments([src], out_dir=tmp_path / "out", runner=parser)

    assert len(parser.calls) == 1
    assert "body from the document" in attachments[0].markdown


def test_docx_is_routed_through_firecrawl_parse(tmp_path):
    src = tmp_path / "spec.docx"
    src.write_bytes(b"PK fake")
    parser = FakeParser()
    load_attachments([src], out_dir=tmp_path / "out", runner=parser)
    assert len(parser.calls) == 1


def test_required_formats_are_all_supported():
    assert ".md" in TEXT_SUFFIXES
    assert {".pdf", ".docx"} <= PARSE_SUFFIXES


# --- parse command --------------------------------------------------------


def test_parse_command_requests_markdown_to_a_file(tmp_path):
    cmd = parse_command(tmp_path / "a.pdf", tmp_path / "out" / "A1.md")
    assert cmd[0] == "firecrawl"
    assert cmd[1] == "parse"
    assert "--format" in cmd and "markdown" in cmd
    assert "-o" in cmd


def test_parse_failure_is_reported_with_the_filename(tmp_path):
    src = tmp_path / "broken.pdf"
    src.write_bytes(b"x")

    def failing(cmd, timeout):
        return 1, "", "unsupported encoding"

    with pytest.raises(AttachmentError, match="broken.pdf"):
        load_attachments([src], out_dir=tmp_path / "out", runner=failing)


def test_empty_parse_output_is_an_error(tmp_path):
    """A parse that yields nothing would silently give members an empty attachment."""
    src = tmp_path / "blank.pdf"
    src.write_bytes(b"x")

    with pytest.raises(AttachmentError):
        load_attachments([src], out_dir=tmp_path / "out", runner=FakeParser(output="   \n"))


# --- labelling and persistence --------------------------------------------


def test_attachments_are_labelled_a1_upwards(tmp_path):
    paths = []
    for name in ("one.md", "two.md", "three.md"):
        p = tmp_path / name
        p.write_text(f"content of {name}")
        paths.append(p)

    attachments = load_attachments(paths, out_dir=tmp_path / "out", runner=FakeParser())
    assert [a.label for a in attachments] == ["A1", "A2", "A3"]


def test_converted_text_is_written_to_the_run_directory(tmp_path):
    src = tmp_path / "notes.md"
    src.write_text("keep me")
    out = tmp_path / "out"
    load_attachments([src], out_dir=out, runner=FakeParser())
    written = list(out.glob("A1*.md"))
    assert written and "keep me" in written[0].read_text()


def test_attachment_keeps_its_original_filename_as_title(tmp_path):
    src = tmp_path / "Quarterly Report.pdf"
    src.write_bytes(b"x")
    attachments = load_attachments([src], out_dir=tmp_path / "out", runner=FakeParser())
    assert attachments[0].title == "Quarterly Report.pdf"


# --- prompt assembly ------------------------------------------------------


def att(label, title, markdown):
    return Attachment(label=label, title=title, source_path=Path(f"/x/{title}"), markdown=markdown)


def test_attachment_block_tags_each_document():
    block = build_attachment_block(
        [att("A1", "one.pdf", "first body"), att("A2", "two.docx", "second body")],
        max_chars_each=1000,
    )
    assert "[A1] one.pdf" in block
    assert "[A2] two.docx" in block
    assert block.index("[A1]") < block.index("[A2]")


def test_attachment_block_truncates_a_long_document():
    block = build_attachment_block([att("A1", "big.pdf", "z" * 50_000)], max_chars_each=500)
    assert "[truncated]" in block
    assert block.count("z") <= 500


def test_attachment_block_is_empty_when_there_are_none():
    assert build_attachment_block([], max_chars_each=1000) == ""
