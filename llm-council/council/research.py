# ABOUTME: Shared web-retrieval stage for the LLM Council.
# ABOUTME: Runs one Firecrawl search, distils the results into a cited brief, and writes raw sources to disk.

"""Stage 0 of the council pipeline.

Every council member reasons over an identical evidence brief so that
disagreement between members reflects interpretation rather than which
documents each happened to retrieve. Raw scrapes are kept alongside the
brief so a member can consult a full source when the brief looks thin.
"""

from __future__ import annotations

import json
import re
import subprocess
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path
from typing import Callable, Iterable, Sequence

from council.prompts import headless_instructions

# agy ignores stdin, so prompts travel as a single argv entry. Linux caps one
# argument at MAX_ARG_STRLEN, which is 32 * PAGE_SIZE and therefore 131072 on
# x86-64. Measured directly here: 131071 bytes execs, 131072 raises E2BIG.
#
# The limit is per argument, not per command line: --model and the rest cost
# nothing against it, and the command line as a whole is bounded separately by
# ARG_MAX (2 MB here), which nothing in this project approaches.
#
# The margin below therefore is not room for other arguments. It exists so an
# oversized prompt fails as a clean CorpusTooLargeError or PromptTooLargeError,
# rather than as an opaque exec failure inside a pane script.
MAX_ARG_STRLEN = 131_072
AGY_ARGV_BUDGET = MAX_ARG_STRLEN - 3_072  # 128_000

# The distiller and the members are bounded by the same kernel limit for the
# same reason, so they share one budget: raising either alone only moves the
# wall from one stage to the next.
DEFAULT_MAX_PROMPT_BYTES = AGY_ARGV_BUDGET
DEFAULT_MAX_CHARS_PER_DOC = 12_000
TRUNCATION_MARKER = "\n\n[truncated]"

# Retrieval aims for this many sources that actually scraped, retrying at
# increasing depth. Firecrawl's search caps at 100 results per call.
DEFAULT_TARGET_SOURCES = 5
DEFAULT_MAX_ATTEMPTS = 15
MAX_SEARCH_LIMIT = 100

# Retrieval is recency-filtered by default: the members are themselves trained on
# historical content, so the marginal value of a live search is what postdates
# their cutoffs. Pass any of NO_RECENCY_VALUES to search all time instead.
DEFAULT_RECENCY_DAYS = 180
NO_RECENCY_VALUES = frozenset({"none", "off", "any", ""})

# agy requires an explicit effort with a bare model name: `--model
# gemini-3.8-flash` alone exits 1. Distillation is mechanical compression
# rather than reasoning, so it runs at the cheapest level.
DISTILLER_EFFORT = "low"

# Stable location the raw scrapes are mirrored to for the escape hatch. Fixed
# filenames because Antigravity matches only literal paths, never globs.
ESCAPE_HATCH_DIR = Path("/config/.llm-council/current/raw")
ESCAPE_HATCH_SLOTS = 10

# A runner takes an argv list and a timeout, and returns (returncode, stdout, stderr).
Runner = Callable[[Sequence[str], int], "tuple[int, str, str]"]


class SearchFailedError(RuntimeError):
    """Firecrawl returned no usable evidence."""


class CorpusTooLargeError(RuntimeError):
    """The assembled prompt would exceed the single-argument limit."""


@dataclass(frozen=True)
class SourceDoc:
    url: str
    title: str
    description: str
    markdown: str
    position: int


@dataclass(frozen=True)
class GatherResult:
    docs: list[SourceDoc]
    credits_used: int
    attempts: int
    target: int


@dataclass(frozen=True)
class ResearchResult:
    brief: str
    docs: list[SourceDoc]
    credits_used: int
    raw_paths: list[Path]
    attempts: int = 1


def _subprocess_runner(cmd: Sequence[str], timeout: int) -> tuple[int, str, str]:
    proc = subprocess.run(
        list(cmd),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    return proc.returncode, proc.stdout, proc.stderr


# --- parsing --------------------------------------------------------------


def parse_search_payload(payload: dict) -> tuple[list[SourceDoc], int]:
    """Turn a `firecrawl search --json` payload into usable source documents.

    Results that carry no markdown failed to scrape and are dropped: an
    untitled URL with no body is not evidence, and passing it to the
    distiller would invite it to speculate about the contents.
    """
    if not payload.get("success"):
        raise SearchFailedError("Firecrawl reported success=false")

    results = (payload.get("data") or {}).get("web") or []
    docs: list[SourceDoc] = []
    for item in results:
        markdown = (item.get("markdown") or "").strip()
        if not markdown:
            continue
        docs.append(
            SourceDoc(
                url=item.get("url", ""),
                title=(item.get("title") or item.get("url") or "Untitled").strip(),
                description=(item.get("description") or "").strip(),
                markdown=markdown,
                position=item.get("position", len(docs) + 1),
            )
        )

    if not docs:
        raise SearchFailedError("Firecrawl returned no scrapeable results")
    return docs, payload.get("creditsUsed", 0)


# --- text helpers ---------------------------------------------------------


def gather_sources(
    query: str,
    target: int = DEFAULT_TARGET_SOURCES,
    max_attempts: int = DEFAULT_MAX_ATTEMPTS,
    tbs: str | None = None,
    out_dir: Path | str = ".",
    runner: Runner = _subprocess_runner,
    timeout: int = 300,
) -> GatherResult:
    """Search repeatedly until `target` sources actually scrape, or attempts run out.

    Firecrawl returns a deterministic ranked list, so re-running the same query at
    the same depth would return the same results. Each retry therefore requests a
    deeper slice and keeps only URLs not already seen.

    Hosts that block scraping (Reddit and LinkedIn among them) return an empty
    body and are skipped, which is the usual reason a nominally successful search
    yields fewer usable sources than requested.
    """
    seen: set[str] = set()
    docs: list[SourceDoc] = []
    credits = 0
    attempts = 0
    limit = target
    out_dir = Path(out_dir)

    while attempts < max_attempts and len(docs) < target:
        attempts += 1
        batch, batch_credits = run_search(
            query,
            limit=limit,
            tbs=tbs,
            out_path=out_dir / f"search-{attempts}.json",
            runner=runner,
            timeout=timeout,
            strict=False,
        )
        credits += batch_credits

        fresh = [doc for doc in batch if doc.url not in seen]
        seen.update(doc.url for doc in batch)
        docs.extend(fresh)

        if len(docs) >= target:
            break

        # At maximum depth with nothing new, further attempts only cost credits.
        if limit >= MAX_SEARCH_LIMIT and not fresh:
            break
        limit = min(MAX_SEARCH_LIMIT, limit * 2)

    if not docs:
        raise SearchFailedError(
            f"no scrapeable sources for {query!r} after {attempts} attempt(s)"
        )

    return GatherResult(
        docs=docs[:target], credits_used=credits, attempts=attempts, target=target
    )


def recency_tbs(today: date, days: int = DEFAULT_RECENCY_DAYS) -> str:
    """Build a Firecrawl `tbs` custom date range covering the last `days` days.

    Uses an explicit range rather than `qdr:m6` so the window is exactly the
    number of days asked for. Dates are unpadded M/D/YYYY, the format the
    Google-style backend behind Firecrawl expects.
    """
    def stamp(d: date) -> str:
        return f"{d.month}/{d.day}/{d.year}"

    start = today - timedelta(days=days)
    return f"cdr:1,cd_min:{stamp(start)},cd_max:{stamp(today)}"


def resolve_tbs(
    tbs: str | None,
    today: date,
    days: int = DEFAULT_RECENCY_DAYS,
) -> str | None:
    """Decide the effective recency filter.

    Recency is the default because the members are themselves trained on
    historical content: the marginal value of retrieval is material published
    after their training cutoffs, and unfiltered search mostly returns older,
    better-linked pages the models already know.

    Note that Firecrawl silently ignores an unrecognised `tbs` rather than
    erroring, so a typo degrades to an unfiltered search without warning.
    """
    if tbs is None:
        return recency_tbs(today, days)
    if tbs.strip().lower() in NO_RECENCY_VALUES:
        return None
    return tbs


def mirror_raw_sources(
    raw_paths: Sequence[Path],
    mirror_dir: Path,
    max_slots: int = ESCAPE_HATCH_SLOTS,
) -> list[Path]:
    """Copy this run's raw scrapes to fixed, stable filenames.

    Antigravity honours only literal path rules in ``permissions.allow`` — no
    glob form matches — so the escape hatch can only work if the readable paths
    never change between runs. Sources land in numbered slots that a fixed set
    of allow-rules can cover.

    Unused slots are overwritten rather than left alone: a shorter run would
    otherwise leave the previous run's sources sitting at a readable path, where
    a member could cite them as current evidence.
    """
    mirror_dir = Path(mirror_dir)
    mirror_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for slot in range(1, max_slots + 1):
        target = mirror_dir / f"S{slot}.md"
        if slot <= len(raw_paths):
            target.write_text(Path(raw_paths[slot - 1]).read_text())
            written.append(target)
        else:
            target.write_text("(unused slot for this run)\n")
    return written


def source_shortfall_warning(requested: int, obtained: int) -> str | None:
    """Warn when fewer sources survived scraping than were asked for.

    Some hosts (Reddit and LinkedIn among them) return an empty body, so a
    nominally successful search can yield far less evidence than requested. A
    council reasoning over one source is worth flagging rather than hiding.
    """
    if obtained >= requested:
        return None
    message = f"only {obtained} of {requested} sources could be scraped"
    if obtained <= 1:
        message += "; the evidence brief will be thin, consider raising --limit"
    return message


def slugify(text: str, maxlen: int = 60) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
    if len(slug) > maxlen:
        slug = slug[:maxlen].rstrip("-")
    return slug or "source"


def truncate_markdown(text: str, max_chars: int) -> tuple[str, bool]:
    """Cap `text`, preferring a line boundary so a source never ends mid-sentence."""
    if len(text) <= max_chars:
        return text, False

    cut = text.rfind("\n", 0, max_chars)
    # Only honour the boundary when it retains most of the budget, otherwise a
    # single long paragraph would collapse the excerpt to almost nothing.
    if cut < max_chars // 2:
        cut = max_chars
    return text[:cut].rstrip() + TRUNCATION_MARKER, True


# --- assembly -------------------------------------------------------------


def build_corpus(docs: Iterable[SourceDoc], max_chars_per_doc: int = DEFAULT_MAX_CHARS_PER_DOC) -> str:
    blocks = []
    for index, source in enumerate(docs, start=1):
        body, _ = truncate_markdown(source.markdown, max_chars_per_doc)
        blocks.append(
            f"[S{index}] {source.title}\nURL: {source.url}\n\n{body}"
        )
    return "\n\n---\n\n".join(blocks)


def build_distill_prompt(
    query: str,
    corpus: str,
    today: date,
    max_prompt_bytes: int = DEFAULT_MAX_PROMPT_BYTES,
) -> str:
    prompt = f"""You are preparing an evidence brief for a panel of three AI models \
who will each answer the same question independently. You are not answering the \
question yourself.

{headless_instructions()}

Today's date is {today.isoformat()}.

QUESTION THE PANEL WILL ANSWER:
{query}

Write a factual brief from the sources below. Requirements:
- Cite every claim with its source tag, for example [S1] or [S2].
- Preserve specifics: version numbers, dates, figures, names, exact quotes.
- Where sources disagree, say so explicitly and cite both.
- Note when a source is undated or looks stale relative to today's date.
- Do not answer the question, recommend anything, or draw conclusions.
- Do not invent sources or cite any tag that is not present below.
- Omit navigation text, marketing copy and boilerplate.
- Start directly with your first section heading. Do not open with a title or a
  date line: both are added around your output.

SOURCES:
{corpus}
"""
    size = len(prompt.encode("utf-8"))
    if size > max_prompt_bytes:
        raise CorpusTooLargeError(
            f"distiller prompt is {size} bytes, over the {max_prompt_bytes} byte limit; "
            "lower --limit or the per-document cap"
        )
    return prompt


def render_raw_document(source: SourceDoc, index: int) -> str:
    """Full, untruncated scrape with front matter, for a member that wants depth."""
    title = source.title.replace('"', "'")
    return (
        "---\n"
        f"source_id: S{index}\n"
        f'title: "{title}"\n'
        f"url: {source.url}\n"
        "---\n\n"
        f"{source.markdown}\n"
    )


def assemble_brief(query: str, distilled: str, docs: Sequence[SourceDoc], today: date) -> str:
    """Attach the canonical source list.

    The URL list is rendered from the Firecrawl payload rather than from the
    distiller's output, so a member can never be pointed at a hallucinated link.
    """
    lines = [
        f"# Evidence brief — {today.isoformat()}",
        "",
        f"**Question:** {query}",
        "",
        distilled.strip(),
        "",
        "## Sources",
        "",
    ]
    lines.extend(
        f"[S{i}] {source.title} — {source.url}" for i, source in enumerate(docs, start=1)
    )
    return "\n".join(lines) + "\n"


# --- process boundaries ---------------------------------------------------


def run_search(
    query: str,
    limit: int,
    tbs: str | None,
    out_path: Path,
    runner: Runner = _subprocess_runner,
    timeout: int = 300,
    strict: bool = True,
) -> tuple[list[SourceDoc], int]:
    """Run one Firecrawl search.

    With ``strict=False`` a batch in which nothing scraped returns no documents
    instead of raising, so a caller escalating through result depths can keep
    going rather than aborting on a shallow miss.

    Firecrawl writes no output file at all when a search matches nothing: it
    prints "No results found." and exits 0. A zero exit therefore does not imply
    the file exists, and a missing file means an empty result set rather than a
    failure to run.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "firecrawl",
        "search",
        query,
        "--scrape",
        "--json",
        "--limit",
        str(limit),
        "--only-main-content",
        "-o",
        str(out_path),
    ]
    if tbs:
        cmd += ["--tbs", tbs]

    code, _, stderr = runner(cmd, timeout)
    if code != 0:
        raise SearchFailedError(f"firecrawl search failed ({code}): {stderr.strip()}")

    if not out_path.exists():
        if strict:
            raise SearchFailedError(
                f"firecrawl returned no results for {query!r} and wrote no output. "
                "A very long or prose-like query is the usual cause; search backends "
                "cap query length well below a full prompt"
            )
        return [], 0

    payload = json.loads(out_path.read_text())
    if strict:
        return parse_search_payload(payload)
    try:
        return parse_search_payload(payload)
    except SearchFailedError:
        return [], payload.get("creditsUsed", 0)


def distil(
    prompt: str,
    model: str,
    runner: Runner = _subprocess_runner,
    timeout: int = 600,
) -> str:
    """Compress the corpus with a cheap model.

    The prompt goes in argv because `agy -p` silently ignores stdin, and a
    zero exit with empty stdout means a tool permission was auto-denied
    rather than that the model had nothing to say.
    """
    cmd = [
        "agy",
        "-p",
        prompt,
        "--model",
        model,
        "--effort",
        DISTILLER_EFFORT,
        "--print-timeout",
        f"{timeout}s",
    ]
    code, stdout, stderr = runner(cmd, timeout + 60)
    if code != 0:
        raise RuntimeError(f"agy distillation failed ({code}): {stderr.strip()}")
    if not stdout.strip():
        raise RuntimeError(
            "agy distillation produced no output; a tool permission was likely auto-denied"
        )
    return stdout.strip()


# --- orchestration --------------------------------------------------------


def research(
    query: str,
    run_dir: Path,
    limit: int = DEFAULT_TARGET_SOURCES,
    max_attempts: int = DEFAULT_MAX_ATTEMPTS,
    tbs: str | None = None,
    distiller_model: str = "gemini-3.8-flash",
    max_chars_per_doc: int = DEFAULT_MAX_CHARS_PER_DOC,
    today: date | None = None,
    search_query: str | None = None,
    mirror_dir: Path = ESCAPE_HATCH_DIR,
    runner: Runner = _subprocess_runner,
) -> ResearchResult:
    """Retrieve, distil and assemble the shared evidence brief.

    ``search_query`` exists because the question a panel is asked and the string
    a search engine can act on are not always the same. A long structured prompt
    makes an unsearchable query, so retrieval can be given its own short one
    while the distiller still sees the real question and keeps what matters.
    """
    today = today or date.today()
    research_dir = run_dir / "research"
    raw_dir = research_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    gathered = gather_sources(
        search_query or query,
        target=limit,
        max_attempts=max_attempts,
        tbs=resolve_tbs(tbs, today),
        out_dir=research_dir,
        runner=runner,
    )
    docs, credits = gathered.docs, gathered.credits_used

    raw_paths = []
    for index, source in enumerate(docs, start=1):
        path = raw_dir / f"S{index}-{slugify(source.title, maxlen=48)}.md"
        path.write_text(render_raw_document(source, index))
        raw_paths.append(path)

    mirror_raw_sources(raw_paths, mirror_dir)

    corpus = build_corpus(docs, max_chars_per_doc=max_chars_per_doc)
    prompt = build_distill_prompt(query, corpus, today)
    distilled = distil(prompt, model=distiller_model, runner=runner)
    brief = assemble_brief(query, distilled, docs, today)

    (research_dir / "brief.md").write_text(brief)
    return ResearchResult(
        brief=brief,
        docs=docs,
        credits_used=credits,
        raw_paths=raw_paths,
        attempts=gathered.attempts,
    )


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="council-research",
        description="Run the council's shared retrieval stage and write an evidence brief.",
    )
    parser.add_argument("query")
    parser.add_argument(
        "--limit",
        type=int,
        default=DEFAULT_TARGET_SOURCES,
        help="target number of sources that actually scrape",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=DEFAULT_MAX_ATTEMPTS,
        help="maximum searches to reach the target",
    )
    parser.add_argument(
        "--tbs",
        default=None,
        help=f"recency filter; defaults to the last {DEFAULT_RECENCY_DAYS} days. "
        "Pass a value such as qdr:w, or 'none' to search all time",
    )
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--distiller-model", default="gemini-3.8-flash")
    args = parser.parse_args(argv)

    result = research(
        args.query,
        run_dir=args.run_dir,
        limit=args.limit,
        max_attempts=args.max_attempts,
        tbs=args.tbs,
        distiller_model=args.distiller_model,
    )
    warning = source_shortfall_warning(args.limit, len(result.docs))
    print(
        f"sources={len(result.docs)} attempts={result.attempts} "
        f"credits={result.credits_used} brief={args.run_dir / 'research' / 'brief.md'}"
    )
    if warning:
        print(f"warning: {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
