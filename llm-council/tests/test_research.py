# ABOUTME: Tests for the shared-retrieval stage of the LLM Council.
# ABOUTME: Covers Firecrawl payload parsing, corpus assembly, argv size guards and brief rendering.

import json
import sys
from datetime import date
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.research import (  # noqa: E402
    DEFAULT_RECENCY_DAYS,
    DISTILLER_EFFORT,
    CorpusTooLargeError,
    SearchFailedError,
    SourceDoc,
    assemble_brief,
    build_corpus,
    build_distill_prompt,
    distil,
    gather_sources,
    mirror_raw_sources,
    parse_search_payload,
    recency_tbs,
    render_raw_document,
    research,
    resolve_tbs,
    run_search,
    slugify,
    source_shortfall_warning,
    truncate_markdown,
)


def doc(n=1, url=None, title=None, markdown="body text"):
    return SourceDoc(
        url=url or f"https://example.com/page-{n}",
        title=title or f"Title {n}",
        description=f"Description {n}",
        markdown=markdown,
        position=n,
    )


# --- parse_search_payload -------------------------------------------------


def test_parse_search_payload_extracts_docs_and_credits():
    payload = {
        "success": True,
        "creditsUsed": 4,
        "data": {
            "web": [
                {
                    "url": "https://a.test/x",
                    "title": "A",
                    "description": "desc a",
                    "markdown": "# A\ncontent",
                    "position": 1,
                },
                {
                    "url": "https://b.test/y",
                    "title": "B",
                    "description": "desc b",
                    "markdown": "# B\ncontent",
                    "position": 2,
                },
            ]
        },
    }
    docs, credits = parse_search_payload(payload)
    assert credits == 4
    assert [d.url for d in docs] == ["https://a.test/x", "https://b.test/y"]
    assert docs[0].title == "A"
    assert docs[0].markdown.startswith("# A")


def test_parse_search_payload_raises_when_unsuccessful():
    with pytest.raises(SearchFailedError):
        parse_search_payload({"success": False, "data": {}})


def test_parse_search_payload_raises_when_no_web_results():
    with pytest.raises(SearchFailedError):
        parse_search_payload({"success": True, "data": {"web": []}})


def test_parse_search_payload_skips_results_without_markdown():
    """A result that failed to scrape carries no markdown and is useless as evidence."""
    payload = {
        "success": True,
        "creditsUsed": 2,
        "data": {
            "web": [
                {"url": "https://a.test", "title": "A", "markdown": "real content"},
                {"url": "https://b.test", "title": "B", "markdown": ""},
                {"url": "https://c.test", "title": "C"},
            ]
        },
    }
    docs, _ = parse_search_payload(payload)
    assert [d.url for d in docs] == ["https://a.test"]


# --- gather_sources (escalating retrieval) --------------------------------


class FakeSearch:
    """Simulates Firecrawl: a fixed ranked result list, some of which fail to scrape.

    Requesting a deeper `--limit` is the only way to surface new URLs, which is
    why retrying with the same limit would be pointless.
    """

    def __init__(self, scrapeable_at, total=100):
        self.scrapeable_at = set(scrapeable_at)
        self.total = total
        self.limits = []

    def __call__(self, cmd, timeout):
        limit = int(cmd[cmd.index("--limit") + 1])
        out_path = Path(cmd[cmd.index("-o") + 1])
        self.limits.append(limit)
        web = [
            {
                "url": f"https://site-{i}.test",
                "title": f"T{i}",
                "markdown": "body" if i in self.scrapeable_at else "",
            }
            for i in range(1, min(limit, self.total) + 1)
        ]
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps({"success": True, "creditsUsed": 2, "data": {"web": web}}))
        return 0, "", ""


def test_gather_stops_after_one_attempt_when_target_is_met(tmp_path):
    fake = FakeSearch(scrapeable_at=range(1, 20))
    result = gather_sources("q", target=5, max_attempts=15, out_dir=tmp_path, runner=fake)
    assert len(result.docs) == 5
    assert result.attempts == 1
    assert fake.limits == [5]


def test_gather_escalates_the_limit_to_reach_deeper_results(tmp_path):
    # Only results 9 and beyond scrape, so limit=5 yields nothing usable.
    fake = FakeSearch(scrapeable_at=range(9, 40))
    result = gather_sources("q", target=5, max_attempts=15, out_dir=tmp_path, runner=fake)
    assert len(result.docs) == 5
    assert fake.limits[0] == 5
    assert fake.limits[1] > fake.limits[0], "a retry must go deeper or it returns the same results"


def test_gather_never_returns_a_duplicate_url(tmp_path):
    fake = FakeSearch(scrapeable_at=range(9, 40))
    result = gather_sources("q", target=5, max_attempts=15, out_dir=tmp_path, runner=fake)
    urls = [d.url for d in result.docs]
    assert len(urls) == len(set(urls))


def test_gather_respects_the_attempt_cap_and_returns_what_it_found(tmp_path):
    fake = FakeSearch(scrapeable_at={1, 2})
    result = gather_sources("q", target=5, max_attempts=3, out_dir=tmp_path, runner=fake)
    assert result.attempts <= 3
    assert len(result.docs) == 2


def test_gather_gives_up_early_once_deeper_searching_cannot_help(tmp_path):
    """At the maximum limit with no new URLs, further attempts only burn credits."""
    fake = FakeSearch(scrapeable_at={1}, total=10)
    result = gather_sources("q", target=5, max_attempts=15, out_dir=tmp_path, runner=fake)
    assert result.attempts < 15


def test_gather_accumulates_credits_across_attempts(tmp_path):
    fake = FakeSearch(scrapeable_at=range(9, 40))
    result = gather_sources("q", target=5, max_attempts=15, out_dir=tmp_path, runner=fake)
    assert result.credits_used == 2 * result.attempts


def test_gather_caps_results_at_the_target(tmp_path):
    """Extra sources would inflate the corpus against agy's argv limit."""
    fake = FakeSearch(scrapeable_at=range(1, 100))
    result = gather_sources("q", target=5, max_attempts=15, out_dir=tmp_path, runner=fake)
    assert len(result.docs) == 5


def test_gather_raises_when_nothing_scrapes_at_all(tmp_path):
    fake = FakeSearch(scrapeable_at=set(), total=10)
    with pytest.raises(SearchFailedError):
        gather_sources("q", target=5, max_attempts=3, out_dir=tmp_path, runner=fake)


# --- recency window -------------------------------------------------------


def test_recency_tbs_spans_exactly_the_requested_window():
    tbs = recency_tbs(date(2026, 7, 27), days=180)
    assert "cd_min:1/28/2026" in tbs
    assert "cd_max:7/27/2026" in tbs
    assert tbs.startswith("cdr:1,")


def test_recency_tbs_uses_unpadded_us_dates():
    """Firecrawl passes tbs through to a Google-style backend: M/D/YYYY, no zero padding."""
    tbs = recency_tbs(date(2026, 3, 5), days=1)
    assert "cd_min:3/4/2026" in tbs
    assert "03/04" not in tbs


def test_default_recency_is_180_days():
    assert DEFAULT_RECENCY_DAYS == 180


def test_resolve_tbs_applies_the_default_window_when_unset():
    resolved = resolve_tbs(None, today=date(2026, 7, 27))
    assert "cd_min:1/28/2026" in resolved


@pytest.mark.parametrize("value", ["none", "None", "off", "any", "", "  "])
def test_resolve_tbs_can_be_switched_off(value):
    """Models are trained on historical content, so recency is the default -- but
    some questions genuinely want the whole archive."""
    assert resolve_tbs(value, today=date(2026, 7, 27)) is None


def test_resolve_tbs_passes_an_explicit_filter_through_untouched():
    assert resolve_tbs("qdr:w", today=date(2026, 7, 27)) == "qdr:w"


# --- stable mirror --------------------------------------------------------


def test_mirror_uses_fixed_slot_filenames(tmp_path):
    """agy only honours literal path rules, so the readable paths must never vary."""
    src = tmp_path / "raw"
    src.mkdir()
    originals = []
    for i in (1, 2):
        p = src / f"S{i}-some-slug.md"
        p.write_text(f"content {i}")
        originals.append(p)

    mirror = tmp_path / "current" / "raw"
    written = mirror_raw_sources(originals, mirror, max_slots=5)

    assert [p.name for p in written] == ["S1.md", "S2.md"]
    assert (mirror / "S1.md").read_text() == "content 1"


def test_mirror_clears_slots_left_over_from_a_previous_run(tmp_path):
    """A shorter run must not leave the last run's sources readable as current."""
    mirror = tmp_path / "current" / "raw"
    mirror.mkdir(parents=True)
    (mirror / "S4.md").write_text("stale content from an older run")

    src = tmp_path / "raw"
    src.mkdir()
    (src / "S1-a.md").write_text("fresh")

    mirror_raw_sources([src / "S1-a.md"], mirror, max_slots=5)

    assert "stale content" not in (mirror / "S4.md").read_text()


def test_mirror_stops_at_the_slot_limit(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    originals = []
    for i in range(1, 9):
        p = src / f"S{i}-x.md"
        p.write_text(f"c{i}")
        originals.append(p)

    written = mirror_raw_sources(originals, tmp_path / "m", max_slots=5)
    assert len(written) == 5


# --- source shortfall -----------------------------------------------------


def test_no_warning_when_every_requested_source_scraped():
    assert source_shortfall_warning(requested=5, obtained=5) is None


def test_warning_names_both_counts_when_sources_are_dropped():
    """Sites like Reddit and LinkedIn return empty markdown; the user should know."""
    warning = source_shortfall_warning(requested=5, obtained=3)
    assert "5" in warning and "3" in warning


def test_warning_is_emphatic_when_evidence_is_very_thin():
    warning = source_shortfall_warning(requested=3, obtained=1)
    assert "thin" in warning.lower()


# --- slugify --------------------------------------------------------------


def test_slugify_produces_filesystem_safe_names():
    assert slugify("Herdr: One Terminal, The Whole Herd!") == "herdr-one-terminal-the-whole-herd"


def test_slugify_truncates_and_never_ends_with_separator():
    out = slugify("a very long title that keeps going and going", maxlen=12)
    assert len(out) <= 12
    assert not out.endswith("-")


def test_slugify_falls_back_when_no_usable_characters():
    assert slugify("???") == "source"


# --- truncate_markdown ----------------------------------------------------


def test_truncate_markdown_leaves_short_text_untouched():
    text, was_truncated = truncate_markdown("short", max_chars=100)
    assert text == "short"
    assert was_truncated is False


def test_truncate_markdown_marks_and_caps_long_text():
    text, was_truncated = truncate_markdown("x" * 500, max_chars=100)
    assert was_truncated is True
    assert "[truncated]" in text
    assert len(text) <= 100 + len("\n\n[truncated]")


def test_truncate_markdown_prefers_a_line_boundary():
    body = "line one\nline two\n" + "z" * 200
    text, _ = truncate_markdown(body, max_chars=30)
    assert "[truncated]" in text
    assert "line one" in text


# --- build_corpus ---------------------------------------------------------


def test_build_corpus_numbers_sources_for_citation():
    corpus = build_corpus([doc(1), doc(2)], max_chars_per_doc=1000)
    assert "[S1]" in corpus
    assert "[S2]" in corpus
    assert "https://example.com/page-1" in corpus
    assert corpus.index("[S1]") < corpus.index("[S2]")


def test_build_corpus_applies_per_document_cap():
    corpus = build_corpus([doc(1, markdown="y" * 5000)], max_chars_per_doc=200)
    assert "[truncated]" in corpus
    assert corpus.count("y") <= 200


# --- build_distill_prompt -------------------------------------------------


def test_distill_prompt_carries_query_date_and_citation_rule():
    prompt = build_distill_prompt(
        query="what changed in X",
        corpus="[S1] ...",
        today=date(2026, 7, 27),
    )
    assert "what changed in X" in prompt
    assert "2026-07-27" in prompt
    assert "[S1]" in prompt
    # The distiller must not invent sources beyond those supplied.
    assert "do not" in prompt.lower()


def test_distill_prompt_forbids_a_duplicate_title_header():
    """assemble_brief supplies the title and date; the distiller must not repeat them."""
    prompt = build_distill_prompt(query="q", corpus="[S1] x", today=date(2026, 7, 27))
    assert "title" in prompt.lower()
    assert "heading" in prompt.lower() or "header" in prompt.lower()


def test_distill_prompt_rejects_corpus_too_large_for_argv():
    with pytest.raises(CorpusTooLargeError):
        build_distill_prompt(
            query="q",
            corpus="z" * 200_000,
            today=date(2026, 7, 27),
            max_prompt_bytes=131_072,
        )


# --- render_raw_document --------------------------------------------------


def test_render_raw_document_has_front_matter_and_full_body():
    out = render_raw_document(doc(3, markdown="full body here"), index=3)
    assert out.startswith("---\n")
    assert "source_id: S3" in out
    assert "url: https://example.com/page-3" in out
    assert "full body here" in out


# --- assemble_brief -------------------------------------------------------


def test_assemble_brief_appends_canonical_sources_not_model_text():
    """URLs come from the payload, never from the distiller, so they cannot be hallucinated."""
    brief = assemble_brief(
        query="q",
        distilled="Some findings citing [S1].",
        docs=[doc(1), doc(2)],
        today=date(2026, 7, 27),
    )
    assert "Some findings citing [S1]." in brief
    assert "2026-07-27" in brief
    assert "[S1] Title 1 — https://example.com/page-1" in brief
    assert "[S2] Title 2 — https://example.com/page-2" in brief


# --- run_search (subprocess boundary, injected runner) --------------------


def test_run_search_builds_expected_firecrawl_command(tmp_path):
    captured = {}

    def fake_runner(cmd, timeout):
        captured["cmd"] = cmd
        out = tmp_path / "out.json"
        payload = {
            "success": True,
            "creditsUsed": 3,
            "data": {"web": [{"url": "https://a.test", "title": "A", "markdown": "body"}]},
        }
        out.write_text(json.dumps(payload))
        return 0, "", ""

    docs, credits = run_search(
        "my query",
        limit=5,
        tbs="qdr:w",
        out_path=tmp_path / "out.json",
        runner=fake_runner,
    )
    cmd = captured["cmd"]
    assert cmd[0] == "firecrawl"
    assert "search" in cmd
    assert "my query" in cmd
    assert "--scrape" in cmd
    assert "--json" in cmd
    assert "--limit" in cmd and "5" in cmd
    assert "--tbs" in cmd and "qdr:w" in cmd
    assert credits == 3
    assert len(docs) == 1


def test_run_search_omits_tbs_when_not_requested(tmp_path):
    captured = {}

    def fake_runner(cmd, timeout):
        captured["cmd"] = cmd
        (tmp_path / "out.json").write_text(
            json.dumps(
                {
                    "success": True,
                    "creditsUsed": 1,
                    "data": {"web": [{"url": "u", "title": "t", "markdown": "m"}]},
                }
            )
        )
        return 0, "", ""

    run_search("q", limit=3, tbs=None, out_path=tmp_path / "out.json", runner=fake_runner)
    assert "--tbs" not in captured["cmd"]


def test_run_search_raises_on_nonzero_exit(tmp_path):
    def fake_runner(cmd, timeout):
        return 1, "", "boom"

    with pytest.raises(SearchFailedError):
        run_search("q", limit=3, tbs=None, out_path=tmp_path / "out.json", runner=fake_runner)


def _no_results_runner(cmd, timeout):
    """Firecrawl on an empty result set: exits 0, says so on stdout, writes no file."""
    return 0, "No results found.\n", ""


def test_run_search_raises_when_firecrawl_wrote_no_file(tmp_path):
    with pytest.raises(SearchFailedError) as exc:
        run_search("q", limit=3, tbs=None, out_path=tmp_path / "out.json", runner=_no_results_runner)
    assert "no results" in str(exc.value).lower()


def test_run_search_treats_a_missing_file_as_empty_when_not_strict(tmp_path):
    docs, credits = run_search(
        "q",
        limit=3,
        tbs=None,
        out_path=tmp_path / "out.json",
        runner=_no_results_runner,
        strict=False,
    )
    assert docs == []
    assert credits == 0


def test_gather_reports_the_query_when_no_search_ever_returned_a_file(tmp_path):
    """An unsearchable query must name itself, not surface a missing-path error."""
    with pytest.raises(SearchFailedError) as exc:
        gather_sources(
            "a query far too long to search",
            target=5,
            out_dir=tmp_path,
            runner=_no_results_runner,
        )
    assert "a query far too long to search" in str(exc.value)


# --- search query separate from the question ------------------------------


class FakeCouncilTools:
    """Stands in for both subprocess boundaries: firecrawl search and agy distil."""

    def __init__(self):
        self.search_queries: list[str] = []
        self.distil_prompts: list[str] = []

    def __call__(self, cmd, timeout):
        if cmd[0] == "firecrawl":
            self.search_queries.append(cmd[2])
            out_path = Path(cmd[cmd.index("-o") + 1])
            out_path.write_text(
                json.dumps(
                    {
                        "success": True,
                        "creditsUsed": 3,
                        "data": {
                            "web": [
                                {"url": "https://a.test", "title": "A", "markdown": "body"}
                            ]
                        },
                    }
                )
            )
            return 0, "", ""
        self.distil_prompts.append(cmd[2])
        return 0, "distilled brief text", ""


def test_research_searches_on_the_search_query_not_the_question(tmp_path):
    """A long prompt is unsearchable, so retrieval must be able to use its own query."""
    tools = FakeCouncilTools()
    research(
        "a very long structured prompt " * 200,
        run_dir=tmp_path,
        limit=1,
        search_query="pptx branding automation",
        today=date(2026, 8, 10),
        mirror_dir=tmp_path / "mirror",
        runner=tools,
    )
    assert tools.search_queries == ["pptx branding automation"]


def test_research_still_shows_the_distiller_the_real_question(tmp_path):
    """The distiller must know what the panel was actually asked, not just the keywords."""
    tools = FakeCouncilTools()
    research(
        "How should I preserve slide masters?",
        run_dir=tmp_path,
        limit=1,
        search_query="python-pptx slide master",
        today=date(2026, 8, 10),
        mirror_dir=tmp_path / "mirror",
        runner=tools,
    )
    assert "How should I preserve slide masters?" in tools.distil_prompts[0]


def test_research_defaults_to_searching_on_the_question(tmp_path):
    tools = FakeCouncilTools()
    research(
        "short answerable question",
        run_dir=tmp_path,
        limit=1,
        today=date(2026, 8, 10),
        mirror_dir=tmp_path / "mirror",
        runner=tools,
    )
    assert tools.search_queries == ["short answerable question"]


# --- distil (subprocess boundary, injected runner) ------------------------


def test_distil_invokes_agy_with_argv_prompt_and_model():
    captured = {}

    def fake_runner(cmd, timeout):
        captured["cmd"] = cmd
        return 0, "the brief text", ""

    out = distil("PROMPT TEXT", model="gemini-3.8-flash", runner=fake_runner)
    cmd = captured["cmd"]
    assert cmd[0] == "agy"
    assert "-p" in cmd
    # Prompt must be passed as argv: agy ignores stdin entirely.
    assert "PROMPT TEXT" in cmd
    assert "--model" in cmd and "gemini-3.8-flash" in cmd
    # agy refuses a bare model name: the effort level is mandatory.
    assert "--effort" in cmd and DISTILLER_EFFORT in cmd
    assert out == "the brief text"


def test_distil_raises_when_agy_fails():
    def fake_runner(cmd, timeout):
        return 1, "", "denied"

    with pytest.raises(RuntimeError):
        distil("p", model="m", runner=fake_runner)


def test_distil_raises_on_empty_output():
    """agy exits 0 while producing nothing when a tool permission is auto-denied."""

    def fake_runner(cmd, timeout):
        return 0, "   \n", ""

    with pytest.raises(RuntimeError):
        distil("p", model="m", runner=fake_runner)
