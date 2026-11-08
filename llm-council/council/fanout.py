# ABOUTME: Stage 1 fan-out: dispatches the same question to every council member in parallel.
# ABOUTME: Members run headless in their own pane and signal completion with a sentinel file.

"""Parallel member dispatch.

Each member runs as a headless one-shot command inside its own Herdr pane, so
the run is visible while it happens, but the orchestrator never reads rendered
terminal output. Answers are captured to files and completion is signalled by a
sentinel, which keeps parsing independent of TUI chrome, ANSI codes and wrapping.

Per-vendor quirks this module exists to absorb:

* ``claude`` reads the prompt from stdin and prints a clean answer to stdout.
* ``codex`` reads stdin via ``-`` but prints token accounting to stdout, so the
  answer must be captured with ``-o``.
* ``agy`` ignores stdin entirely, so its prompt is substituted into argv and is
  therefore bounded by ``MAX_ARG_STRLEN``.
"""

from __future__ import annotations

import re
import shlex
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

from council.research import AGY_ARGV_BUDGET, ESCAPE_HATCH_DIR
from council.prompts import headless_instructions
from council.panes import (
    CHAIRMAN,
    HERDR,
    MEMBER_A,
    MEMBER_B,
    MEMBER_C,
    Runner,
    _subprocess_runner,
    herdr_ok,
)

# agy substitutes the prompt into a single argv entry, so it is bounded by the
# kernel's per-argument limit. Shared with the distiller: see AGY_ARGV_BUDGET in
# research.py for why the margin is what it is.
AGY_ARGV_LIMIT = AGY_ARGV_BUDGET

# Effort is set explicitly rather than left to agy's default. Members reason at
# full stretch, matching codex's reasoning effort.
AGY_EFFORT = "high"

# Members that can actually read a file here. Codex is excluded: its filesystem
# sandbox fails closed in this container, so it cannot read even its own working
# directory. Offering it the escape hatch would just waste a turn.
HATCH_CAPABLE_KINDS = frozenset({"claude", "agy"})

# Self-reported marker a member appends when it opens a full source, so hatch
# usage -- and the asymmetry it creates -- stays visible in the run record.
HATCH_MARKER = re.compile(r"\[S(\d+):full\]")


class PromptTooLargeError(RuntimeError):
    """A prompt exceeds the argv limit of an argv-bound agent."""


class StageTimeout(RuntimeError):
    """One or more members did not finish within the allotted time."""


class AgentFailedError(RuntimeError):
    """A member exited non-zero or produced nothing usable."""


@dataclass(frozen=True)
class Member:
    role: str
    label: str
    kind: str
    model: str


MEMBERS: tuple[Member, ...] = (
    Member(MEMBER_A, "A", "claude", "claude-opus-5-5"),
    Member(MEMBER_B, "B", "codex", "gpt-6-astra"),
    Member(MEMBER_C, "C", "agy", "gemini-4-argon"),
)

MEMBERS_BY_ROLE = {member.role: member for member in MEMBERS}

# The chair reuses the member machinery: same script rendering, same sentinel
# protocol, its own pane.
CHAIRMAN_MEMBER = Member(CHAIRMAN, "chair", "claude", "claude-fable-5-1")

# Optional `--debate` mode. Stance is fixed per pane rather than shuffled so the
# layout stays legible: pane A is always the advocate. This deliberately trades
# the "which model is best" comparison for coverage of the argument space.
DEBATE_STANCES: dict[str, str] = {
    MEMBER_A: (
        "Build the strongest defensible case for the answer you think is right. "
        "Do not hedge it into vagueness; make the best version of the argument."
    ),
    MEMBER_B: (
        "Attack the question's easy answer. Surface failure modes, hidden costs, "
        "and the conditions under which the obvious choice is wrong."
    ),
    MEMBER_C: (
        "Stay even-handed. Weigh the considerations on both sides and say which "
        "way the balance falls and how confidently."
    ),
}


@dataclass(frozen=True)
class MemberPaths:
    prompt: Path
    script: Path
    answer: Path
    status: Path
    done: Path
    log: Path


def member_paths(run_dir: Path | str, role: str, stage: str) -> MemberPaths:
    base = Path(run_dir) / stage
    return MemberPaths(
        prompt=base / f"{role}.prompt.md",
        script=base / f"{role}.sh",
        answer=base / f"{role}.md",
        status=base / f"{role}.status",
        done=base / f"{role}.done",
        log=base / f"{role}.log",
    )


# --- prompts --------------------------------------------------------------


def build_stage1_prompt(
    question: str,
    brief: str | None,
    stance: str | None = None,
    escape_hatch_slots: int = 0,
    attachments: str | None = None,
) -> str:
    sections = [
        "Answer the question below as well as you can.",
        "",
        headless_instructions(allow_source_reads=bool(brief) and escape_hatch_slots > 0),
        "",
        "Write the answer only. No preamble, no restatement of the question, and no "
        "sign-off. Be concrete and state your reasoning. Where you are uncertain, say "
        "so plainly rather than hedging everything.",
        "",
        f"QUESTION:\n{question}",
    ]
    if attachments:
        sections += [
            "",
            "The user attached the documents below. Every member of this panel received "
            "the same attachments. Treat them as primary evidence, ahead of your own "
            "recollection, and cite them by tag, for example [A1]. A document marked "
            "[truncated] has had its tail cut for length; say so rather than guessing "
            "at what it contained.",
            "",
            f"ATTACHED DOCUMENTS:\n{attachments}",
        ]
    if stance:
        sections += ["", f"YOUR ROLE ON THIS PANEL:\n{stance}"]
    if brief:
        sections += [
            "",
            "You have been given a shared evidence brief. Every member of this panel "
            "received exactly the same brief. Ground factual claims in it and cite its "
            "source tags, for example [S1].",
            "",
            "Two rules about the brief versus your own knowledge, and they pull in "
            "opposite directions on purpose:",
            "",
            "1. Where the brief and your own knowledge disagree on a fact, the brief "
            "wins. It was retrieved recently and your training may predate it. Say "
            "explicitly when this happens and what you previously believed — a "
            "contradiction is useful information, not something to smooth over.",
            "",
            "2. Where the brief is silent or thin, use your own knowledge rather than "
            "pretending the question is unanswerable. The brief is a handful of scraped "
            "pages, not a survey, and its coverage may be poor or skewed. Mark any claim "
            "resting on your own knowledge rather than a source with [prior], so the "
            "reader can tell retrieved fact from recollection.",
            "",
            "Do not guess. If you know neither the brief nor your own knowledge settles "
            "something, say that plainly.",
            "",
            f"EVIDENCE BRIEF:\n{brief}",
        ]
        if escape_hatch_slots > 0:
            # Literal paths only: Antigravity's permission matcher rejects globs,
            # and a member that greps for these files is denied.
            paths = "\n".join(
                f"  {ESCAPE_HATCH_DIR}/S{i}.md" for i in range(1, escape_hatch_slots + 1)
            )
            sections += [
                "",
                "If the brief is too compressed to settle something, you may read the "
                "full text of a source. Source [Sn] is at the exact path below. Open a "
                "path directly; do not search or glob for these files. Whenever you use "
                "a claim drawn from a full source rather than from the brief, cite it as "
                "[Sn:full] so the record shows which sources were opened.",
                "",
                paths,
            ]
    return "\n".join(sections) + "\n"


# --- script rendering -----------------------------------------------------


def _agent_command(
    member: Member, paths: MemberPaths, timeout_s: int, escape_hatch: bool = False
) -> str:
    # Every command pipes through `tee` so the pane shows the work as it happens
    # while the orchestrator still reads a clean file. `pipefail` is set, so the
    # agent's exit status survives the pipe.
    if member.kind == "claude":
        # stdout is clean; tools are disabled so the member answers rather than explores.
        # With the escape hatch on, Read alone is granted -- never Edit, Write or Bash.
        tools = "Read" if (escape_hatch and member.kind in HATCH_CAPABLE_KINDS) else ""
        return (
            f"claude -p --model {member.model} --allowedTools \"{tools}\" "
            f"< {shlex.quote(str(paths.prompt))} 2> {shlex.quote(str(paths.log))} "
            f"| tee {shlex.quote(str(paths.answer))}"
        )

    if member.kind == "codex":
        # stdout carries token accounting, so -o is the only clean capture;
        # that frees stdout to be the visible progress stream.
        return (
            f"codex exec - --skip-git-repo-check --sandbox read-only "
            f"--model {member.model} "
            f"-c model_reasoning_effort=high "
            f"-o {shlex.quote(str(paths.answer))} "
            f"< {shlex.quote(str(paths.prompt))} 2>&1 "
            f"| tee {shlex.quote(str(paths.log))}"
        )

    if member.kind == "agy":
        # agy ignores stdin, so the prompt is substituted into argv.
        return (
            f"agy -p \"$(cat {shlex.quote(str(paths.prompt))})\" "
            f"--model {member.model} --effort {AGY_EFFORT} "
            f"--print-timeout {timeout_s}s "
            f"2> {shlex.quote(str(paths.log))} "
            f"| tee {shlex.quote(str(paths.answer))}"
        )

    raise ValueError(f"unsupported member kind: {member.kind}")


def render_member_script(
    member: Member, paths: MemberPaths, timeout_s: int, escape_hatch: bool = False
) -> str:
    """Build the shell script a pane runs for this member.

    Deliberately omits ``set -e``: an aborting script would never reach the
    sentinel, and the orchestrator would wait on it until the stage timed out.
    """
    if member.kind == "agy" and paths.prompt.exists():
        size = len(paths.prompt.read_text().encode("utf-8"))
        if size > AGY_ARGV_LIMIT:
            raise PromptTooLargeError(
                f"{member.label}: prompt is {size} bytes, over the {AGY_ARGV_LIMIT} byte "
                "argv limit for agy; shorten the brief or lower the source count"
            )

    return f"""#!/usr/bin/env bash
# Council member {member.label} ({member.kind} / {member.model}).
# No `set -e`: the sentinel must fire even when the agent fails.
set -uo pipefail

printf '\\n=== member {member.label} · {member.kind} {member.model} ===\\n'

{_agent_command(member, paths, timeout_s, escape_hatch)}
status=$?

echo "$status" > {shlex.quote(str(paths.status))}
printf 'member {member.label} finished with status %s\\n' "$status"
touch {shlex.quote(str(paths.done))}
"""


# --- dispatch -------------------------------------------------------------


def dispatch_stage(
    panes: Mapping[str, str],
    run_dir: Path | str,
    stage: str,
    prompts: Mapping[str, str],
    runner: Runner = _subprocess_runner,
    timeout_s: int = 900,
    members: Sequence[Member] = MEMBERS,
    escape_hatch: bool = False,
) -> None:
    """Write each member's prompt and script, then start them all in their panes."""
    run_dir = Path(run_dir)
    (run_dir / stage).mkdir(parents=True, exist_ok=True)

    for member in members:
        paths = member_paths(run_dir, member.role, stage)
        paths.prompt.write_text(prompts[member.role])

        # A sentinel left over from an earlier run would be read as this run's result.
        for stale in (paths.done, paths.status, paths.answer, paths.log):
            stale.unlink(missing_ok=True)

        paths.script.write_text(
            render_member_script(member, paths, timeout_s, escape_hatch)
        )

    for member in members:
        paths = member_paths(run_dir, member.role, stage)
        herdr_ok(
            [HERDR, "pane", "run", panes[member.role], "bash", str(paths.script)],
            runner=runner,
        )


def wait_for_stage(
    roles: Iterable[str],
    run_dir: Path | str,
    stage: str,
    timeout: float = 900,
    poll: float = 1.0,
    clock: Callable[[], float] = time.monotonic,
    sleep: Callable[[float], None] = time.sleep,
) -> dict[str, int]:
    """Block until every member's sentinel appears, returning their exit statuses."""
    roles = list(roles)
    deadline = clock() + timeout
    statuses: dict[str, int] = {}

    while True:
        for role in roles:
            if role in statuses:
                continue
            paths = member_paths(run_dir, role, stage)
            if paths.done.exists():
                raw = paths.status.read_text().strip() if paths.status.exists() else ""
                statuses[role] = int(raw) if raw.lstrip("-").isdigit() else -1

        if len(statuses) == len(roles):
            return statuses

        if clock() >= deadline:
            stalled = [role for role in roles if role not in statuses]
            raise StageTimeout(
                f"{stage} timed out after {timeout}s; still waiting on: {', '.join(stalled)}"
            )
        sleep(poll)


def collect_answers(roles: Iterable[str], run_dir: Path | str, stage: str) -> dict[str, str]:
    answers: dict[str, str] = {}
    for role in roles:
        paths = member_paths(run_dir, role, stage)
        status = paths.status.read_text().strip() if paths.status.exists() else "missing"
        log_tail = paths.log.read_text().strip()[-500:] if paths.log.exists() else ""

        if status != "0":
            raise AgentFailedError(f"{role} exited with status {status}: {log_tail}")

        text = paths.answer.read_text().strip() if paths.answer.exists() else ""
        if not text:
            raise AgentFailedError(
                f"{role} exited cleanly but produced no answer "
                f"(a tool permission may have been auto-denied): {log_tail}"
            )
        answers[role] = text
    return answers


def detect_hatch_usage(answers: Mapping[str, str]) -> dict[str, list[int]]:
    """Which sources each member reports having opened in full.

    Self-reported, so it undercounts a member that forgets the marker. It exists
    to make the Claude/Gemini-vs-Codex asymmetry visible, not to audit it.
    """
    return {
        role: sorted({int(n) for n in HATCH_MARKER.findall(text or "")})
        for role, text in answers.items()
    }
