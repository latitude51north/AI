# ABOUTME: End-to-end orchestration of a council run across all three stages.
# ABOUTME: Wires retrieval, fan-out, peer review and synthesis into one command.

"""The full pipeline.

    stage 0  one shared Firecrawl retrieval, distilled into a cited brief
    stage 1  three members answer the same question in parallel
    stage 2  each member ranks all three answers, anonymised and shuffled
    stage 3  the chair synthesises a single answer

Everything is written to a run directory so a run can be inspected, diffed or
re-read afterwards. Members communicate only through files; no stage depends on
reading rendered terminal output.
"""

from __future__ import annotations

import json
import random
import re
import sys
from datetime import date, datetime
from pathlib import Path
from typing import Mapping, Sequence

from council.attachments import (
    build_attachment_block,
    load_attachments,
    validate_attachment_paths,
)
from council.deliberate import (
    MalformedRankingError,
    aggregate_rankings,
    assign_labels,
    build_chair_prompt,
    build_review_prompt,
    divergence_report,
    parse_ranking,
)
from council.export import ExportError, export_document, validate_export_target
from council.fanout import (
    CHAIRMAN_MEMBER,
    DEBATE_STANCES,
    HATCH_CAPABLE_KINDS,
    MEMBERS,
    detect_hatch_usage,
    build_stage1_prompt,
    collect_answers,
    dispatch_stage,
    wait_for_stage,
)
from council.panes import (
    CHAIRMAN,
    DEFAULT_TAB_LABEL,
    HERDR,
    ORCHESTRATOR,
    PaneSet,
    Runner,
    _subprocess_runner,
    close_layout,
    create_layout,
    find_tab_by_label,
    herdr_ok,
)
from council.research import ESCAPE_HATCH_DIR
from council import vault

RUNS_ROOT = Path("/config/.llm-council/runs")


# --- helpers --------------------------------------------------------------


def default_run_dir(question: str, now: datetime | None = None, root: Path = RUNS_ROOT) -> Path:
    now = now or datetime.now()
    slug = re.sub(r"[^a-z0-9]+", "-", question.lower()).strip("-")[:48].rstrip("-") or "query"
    return root / f"{now:%Y%m%d-%H%M%S}-{slug}"


def parse_ballots(
    reviews: Mapping[str, str],
    expected: set[str],
    labels: Mapping[str, str],
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Parse each reviewer's ballot, keeping the usable ones.

    Ballots are keyed by the reviewer's own label so nothing downstream can map
    a ballot back to a vendor. A reviewer that produced no parseable ranking is
    recorded as a failure rather than aborting the stage.
    """
    ballots: dict[str, list[str]] = {}
    failures: dict[str, str] = {}
    for role, text in reviews.items():
        try:
            ballots[labels[role]] = parse_ranking(text, expected=expected)
        except MalformedRankingError as exc:
            failures[role] = str(exc)
    return ballots, failures


def build_final_document(
    question: str,
    chair_answer: str,
    order: Sequence[str],
    divergence: Mapping,
    labels: Mapping[str, str],
    members: Mapping[str, str],
) -> str:
    """The chair's answer plus provenance.

    The label-to-vendor mapping is disclosed here: anonymity exists to stop
    members favouring each other, not to keep the reader in the dark.
    """
    lines = [
        f"# {question}",
        "",
        chair_answer.strip(),
        "",
        "---",
        "",
        "## How this answer was produced",
        "",
        "| Shown to panel as | Member |",
        "| --- | --- |",
    ]
    for role, label in sorted(labels.items(), key=lambda kv: kv[1]):
        lines.append(f"| {label} | {members.get(role, role)} |")

    lines += [
        "",
        f"- Panel ranking: {' > '.join(order) if order else 'no usable ballots'}",
        f"- Answer mean similarity: {divergence.get('mean_similarity')}",
    ]
    if divergence.get("converged"):
        lines.append(
            "- **Members converged.** They said much the same thing, so the council added "
            "little over a single call on this question."
        )
    return "\n".join(lines) + "\n"


def _ranked_members(
    order: Sequence[str], labels: Mapping[str, str], members: Mapping[str, str]
) -> list[str]:
    """Translate a ranking of panel letters back into member names."""
    role_by_label = {label: role for role, label in labels.items()}
    return [members.get(role_by_label[label], label) for label in order]


def _announce(pane: str | None, message: str, runner: Runner = _subprocess_runner) -> None:
    """Echo progress into the orchestrator pane so the layout tells the whole story."""
    if not pane:
        return
    try:
        herdr_ok([HERDR, "pane", "run", pane, "echo", message], runner=runner)
    except Exception:  # noqa: BLE001 - progress display must never fail a run
        pass


# --- pipeline -------------------------------------------------------------


def run_council(
    question: str,
    run_dir: Path,
    web: bool = False,
    limit: int = 5,
    max_attempts: int = 15,
    tbs: str | None = None,
    search_query: str | None = None,
    debate: bool = False,
    escape_hatch: bool = False,
    attach: Sequence[Path] = (),
    export: Path | None = None,
    vault_dir: Path | None = None,
    seed: int | None = None,
    timeout: float = 900,
    tab_label: str = DEFAULT_TAB_LABEL,
    keep_panes: bool = True,
    mirror_dir: Path = ESCAPE_HATCH_DIR,
    runner: Runner = _subprocess_runner,
) -> Path:
    """Drive a whole council run and return the path to the final document.

    ``runner`` is the single seam through which every child process is started —
    Herdr, Firecrawl, agy and uv alike. Injecting it lets the whole pipeline be
    exercised, stage sequencing included, without spending credits or quota.
    """
    # Validate attachments before retrieval so a bad path fails before any
    # credits or subscription quota are spent.
    validate_attachment_paths(attach)
    if export:
        validate_export_target(export)

    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "query.txt").write_text(question + "\n")

    roles = [m.role for m in MEMBERS]
    member_names = {m.role: f"{m.kind} {m.model}" for m in MEMBERS}

    attachment_block = None
    if attach:
        loaded = load_attachments(attach, out_dir=run_dir / "attachments", runner=runner)
        attachment_block = build_attachment_block(loaded)
        print(f"attachments: {len(loaded)} document(s) converted", flush=True)

    brief = None
    source_count = 0
    if web:
        from council.research import research, source_shortfall_warning

        effective_query = search_query or question
        shown = effective_query if len(effective_query) <= 70 else effective_query[:67] + "..."
        print(f"stage 0: retrieving on {shown!r}...", flush=True)
        result = research(
            question,
            run_dir=run_dir,
            limit=limit,
            max_attempts=max_attempts,
            tbs=tbs,
            search_query=search_query,
            today=date.today(),
            mirror_dir=mirror_dir,
            runner=runner,
        )
        brief = result.brief
        source_count = len(result.docs)
        print(
            f"stage 0: {len(result.docs)} sources in {result.attempts} attempt(s), "
            f"{result.credits_used} credits",
            flush=True,
        )
        warning = source_shortfall_warning(limit, len(result.docs))
        if warning:
            print(f"stage 0: warning, {warning}", flush=True)

    panes = create_layout(label=tab_label, cwd=run_dir, focus=False, runner=runner)
    panes.save(run_dir / "panes.json")
    orchestrator = panes.roles.get(ORCHESTRATOR)
    _announce(orchestrator, f"council: {question}", runner=runner)

    # Stage 1 -------------------------------------------------------------
    print("stage 1: members answering...", flush=True)
    _announce(orchestrator, "stage 1/3: members answering", runner=runner)
    prompts = {
        m.role: build_stage1_prompt(
            question,
            brief=brief,
            stance=DEBATE_STANCES[m.role] if debate else None,
            escape_hatch_slots=(
                source_count
                if escape_hatch and m.kind in HATCH_CAPABLE_KINDS
                else 0
            ),
            attachments=attachment_block,
        )
        for m in MEMBERS
    }
    dispatch_stage(
        panes=panes.roles,
        run_dir=run_dir,
        stage="stage1",
        prompts=prompts,
        escape_hatch=escape_hatch,
        runner=runner,
    )
    wait_for_stage(roles, run_dir, "stage1", timeout=timeout)
    answers = collect_answers(roles, run_dir, "stage1")

    if escape_hatch:
        usage = detect_hatch_usage(answers)
        (run_dir / "hatch.json").write_text(json.dumps(usage, indent=2))
        opened = {role: srcs for role, srcs in usage.items() if srcs}
        print(f"stage 1: full sources opened by {len(opened)} member(s)", flush=True)

    labels = assign_labels(roles, rng=random.Random(seed))
    (run_dir / "labels.json").write_text(json.dumps(labels, indent=2))

    divergence = divergence_report(answers, labels)
    (run_dir / "divergence.json").write_text(json.dumps(divergence, indent=2))
    print(f"stage 1: mean similarity {divergence['mean_similarity']}", flush=True)

    # Stage 2 -------------------------------------------------------------
    print("stage 2: peer review...", flush=True)
    _announce(orchestrator, "stage 2/3: peer review", runner=runner)
    review_prompt = build_review_prompt(
        question, answers, labels, brief=brief, attachments=attachment_block
    )
    dispatch_stage(
        panes=panes.roles,
        run_dir=run_dir,
        stage="stage2",
        prompts={role: review_prompt for role in roles},
        runner=runner,
    )
    wait_for_stage(roles, run_dir, "stage2", timeout=timeout)
    reviews = collect_answers(roles, run_dir, "stage2")

    ballots, failures = parse_ballots(reviews, expected=set(labels.values()), labels=labels)
    order, scores = aggregate_rankings(ballots)
    if not order:
        order = sorted(labels.values())
    (run_dir / "rankings.json").write_text(
        json.dumps(
            {"ballots": ballots, "scores": scores, "order": order, "failures": failures},
            indent=2,
        )
    )
    if failures:
        print(f"stage 2: {len(failures)} unusable ballot(s)", flush=True)
    print(f"stage 2: panel ranking {' > '.join(order)}", flush=True)

    # Stage 3 -------------------------------------------------------------
    print("stage 3: chair synthesising...", flush=True)
    _announce(orchestrator, "stage 3/3: chair synthesising", runner=runner)
    chair_prompt = build_chair_prompt(
        question, answers, labels, order=order, brief=brief, rankings=ballots
    )
    dispatch_stage(
        panes=panes.roles,
        run_dir=run_dir,
        stage="stage3",
        prompts={CHAIRMAN: chair_prompt},
        members=[CHAIRMAN_MEMBER],
        runner=runner,
    )
    wait_for_stage([CHAIRMAN], run_dir, "stage3", timeout=timeout)
    chair_answer = collect_answers([CHAIRMAN], run_dir, "stage3")[CHAIRMAN]

    final_path = run_dir / "final.md"
    final_path.write_text(
        build_final_document(question, chair_answer, order, divergence, labels, member_names)
    )
    if export:
        # A failed conversion must not discard a council run that succeeded.
        try:
            written = export_document(final_path, export, runner=runner)
            print(f"exported: {written}", flush=True)
        except ExportError as exc:
            print(f"export failed, {final_path} is still intact: {exc}", file=sys.stderr)

    # The vault folder is resolved here rather than as a default argument so
    # tests can redirect it away from the real vault.
    try:
        note = vault.write_vault_note(
            final_path.read_text(),
            question,
            run_dir=run_dir,
            panel=[member_names[role] for role in roles],
            chair=f"{CHAIRMAN_MEMBER.kind} {CHAIRMAN_MEMBER.model}",
            ranking=_ranked_members(order, labels, member_names) if ballots else [],
            mean_similarity=divergence.get("mean_similarity"),
            vault_dir=vault_dir or vault.VAULT_DIR,
        )
        print(f"note: {note}", flush=True)
    except (vault.VaultError, OSError) as exc:
        print(f"vault note failed, {final_path} is still intact: {exc}", file=sys.stderr)

    _announce(orchestrator, f"done: {final_path}", runner=runner)

    if not keep_panes:
        tab_id = find_tab_by_label(tab_label, runner=runner)
        if tab_id:
            close_layout(tab_id, runner=runner)

    return final_path


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="council run",
        description="Ask a question to a panel of three models and synthesise their answers.",
    )
    parser.add_argument("question")
    parser.add_argument("--web", action="store_true", help="retrieve a shared evidence brief first")
    parser.add_argument(
        "--limit", type=int, default=5, help="target sources that actually scrape, with --web"
    )
    parser.add_argument(
        "--max-attempts", type=int, default=15, help="maximum searches to reach --limit"
    )
    parser.add_argument(
        "--tbs",
        default=None,
        help="recency filter for --web; defaults to the last 180 days. "
        "Pass a value such as qdr:w, or 'none' to search all time",
    )
    parser.add_argument(
        "--query",
        default=None,
        metavar="TEXT",
        help="search on this instead of the question, for --web. Use when the question "
        "is a long prompt: search backends cap query length well below a full prompt",
    )
    parser.add_argument(
        "--debate",
        action="store_true",
        help="give each member a different stance: advocate, critic, neutral",
    )
    parser.add_argument(
        "--escape-hatch",
        action="store_true",
        help="let members read full sources when the brief is thin (needs --web)",
    )
    parser.add_argument(
        "--attach",
        type=Path,
        action="append",
        default=[],
        metavar="PATH",
        help="attach a local document (md, txt, pdf, docx, doc, odt, rtf, xlsx, xls, html). "
        "Repeatable",
    )
    parser.add_argument(
        "--export",
        type=Path,
        default=None,
        metavar="PATH",
        help="also write the synthesis here; .md copies it, .docx converts via pandoc",
    )
    parser.add_argument("--seed", type=int, default=None, help="fix the anonymisation shuffle")
    parser.add_argument("--run-dir", type=Path, default=None)
    parser.add_argument("--timeout", type=float, default=900, help="per-stage timeout in seconds")
    parser.add_argument("--close-panes", action="store_true", help="tear the layout down when done")
    args = parser.parse_args(argv)

    # Without --web nothing is retrieved, so a search query would be silently
    # discarded. Fail loudly rather than let it look as though it took effect.
    if args.query and not args.web:
        parser.error("--query only applies with --web")

    run_dir = args.run_dir or default_run_dir(args.question)
    try:
        final_path = run_council(
            args.question,
            run_dir=run_dir,
            web=args.web,
            limit=args.limit,
            max_attempts=args.max_attempts,
            tbs=args.tbs,
            search_query=args.query,
            debate=args.debate,
            escape_hatch=args.escape_hatch,
            attach=args.attach,
            export=args.export,
            seed=args.seed,
            timeout=args.timeout,
            keep_panes=not args.close_panes,
        )
    except Exception as exc:  # noqa: BLE001 - surface the cause, not a traceback
        print(f"council run failed: {exc}", file=sys.stderr)
        return 1

    print(f"\n{final_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
