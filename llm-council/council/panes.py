# ABOUTME: Builds and tears down the Herdr pane layout that hosts the council.
# ABOUTME: Pure layout plan plus a thin wrapper over the Herdr socket API.

"""Stage scaffolding.

Each council role gets its own pane so the run is visible while it happens.
The panes are plain terminals: council stages are dispatched into them later
as headless commands, and answers are captured from files rather than scraped
from the terminal, so nothing here depends on reading rendered output.

Herdr's ``--ratio`` is the fraction of the split axis retained by the *source*
pane; the new pane receives the remainder.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Sequence

HERDR = "herdr"

ORCHESTRATOR = "orchestrator"
MEMBER_A = "member_a"
MEMBER_B = "member_b"
MEMBER_C = "member_c"
CHAIRMAN = "chairman"

ROLES = [ORCHESTRATOR, MEMBER_A, MEMBER_B, MEMBER_C, CHAIRMAN]

PANE_LABELS = {
    ORCHESTRATOR: "orchestrator",
    MEMBER_A: "A · claude opus 5.5",
    MEMBER_B: "B · codex gpt-6-astra",
    MEMBER_C: "C · gemini 4 argon",
    CHAIRMAN: "chair · claude fable 5.1",
}

DEFAULT_TAB_LABEL = "council"

Runner = Callable[[Sequence[str], int], "tuple[int, str, str]"]


class HerdrError(RuntimeError):
    """The Herdr socket API rejected a command or returned something unusable."""


@dataclass(frozen=True)
class SplitStep:
    source_role: str
    new_role: str
    direction: str
    ratio: float


# Orchestrator strip on top, then a 2x2 grid: members A and B above, member C
# and the chairman below. Halves stay readable where thirds did not.
# member_a starts as the whole body and is subdivided by the later steps.
LAYOUT_PLAN: list[SplitStep] = [
    SplitStep(ORCHESTRATOR, MEMBER_A, "down", 0.22),
    SplitStep(MEMBER_A, MEMBER_C, "down", 0.50),
    SplitStep(MEMBER_A, MEMBER_B, "right", 0.50),
    SplitStep(MEMBER_C, CHAIRMAN, "right", 0.50),
]


@dataclass(frozen=True)
class PaneSet:
    tab_id: str
    roles: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        missing = set(ROLES) - set(self.roles)
        if missing:
            raise ValueError(f"pane set is missing roles: {sorted(missing)}")

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"tab_id": self.tab_id, "roles": self.roles}, indent=2))

    @classmethod
    def load(cls, path: Path) -> "PaneSet":
        data = json.loads(Path(path).read_text())
        return cls(tab_id=data["tab_id"], roles=data["roles"])


def _subprocess_runner(cmd: Sequence[str], timeout: int) -> tuple[int, str, str]:
    proc = subprocess.run(list(cmd), capture_output=True, text=True, timeout=timeout)
    return proc.returncode, proc.stdout, proc.stderr


def herdr_json(cmd: Sequence[str], runner: Runner = _subprocess_runner, timeout: int = 30) -> dict:
    """Run a Herdr command and unwrap its ``{"id": ..., "result": ...}`` envelope."""
    code, stdout, stderr = runner(cmd, timeout)
    if code != 0:
        raise HerdrError(f"{' '.join(cmd)} failed ({code}): {stderr.strip()}")
    try:
        payload = json.loads(stdout)
    except json.JSONDecodeError as exc:
        raise HerdrError(f"{' '.join(cmd)} returned non-JSON output: {stdout[:200]!r}") from exc
    if "result" not in payload:
        raise HerdrError(f"{' '.join(cmd)} returned no result: {stdout[:200]!r}")
    return payload["result"]


def herdr_ok(cmd: Sequence[str], runner: Runner = _subprocess_runner, timeout: int = 30) -> None:
    """Run a Herdr command that reports success by exit code alone.

    Not every command returns an envelope: ``pane run`` succeeds silently, so
    demanding JSON from it would turn every successful dispatch into an error.
    """
    code, _, stderr = runner(cmd, timeout)
    if code != 0:
        raise HerdrError(f"{' '.join(cmd)} failed ({code}): {stderr.strip()}")


def _focus_flag(focus: bool) -> str:
    return "--focus" if focus else "--no-focus"


def find_tab_by_label(label: str, runner: Runner = _subprocess_runner) -> str | None:
    result = herdr_json([HERDR, "tab", "list"], runner=runner)
    for tab in result.get("tabs", []):
        if tab.get("label") == label:
            return tab.get("tab_id")
    return None


def create_layout(
    label: str = DEFAULT_TAB_LABEL,
    cwd: str | Path = ".",
    runner: Runner = _subprocess_runner,
    focus: bool = False,
) -> PaneSet:
    """Create the council tab and split it into one pane per role."""
    cwd = str(cwd)
    focus_flag = _focus_flag(focus)

    created = herdr_json(
        [HERDR, "tab", "create", "--label", label, "--cwd", cwd, focus_flag],
        runner=runner,
    )
    tab_id = created["tab"]["tab_id"]
    roles = {ORCHESTRATOR: created["root_pane"]["pane_id"]}

    for step in LAYOUT_PLAN:
        result = herdr_json(
            [
                HERDR,
                "pane",
                "split",
                roles[step.source_role],
                "--direction",
                step.direction,
                "--ratio",
                str(step.ratio),
                "--cwd",
                cwd,
                focus_flag,
            ],
            runner=runner,
        )
        roles[step.new_role] = result["pane"]["pane_id"]

    for role in ROLES:
        herdr_json([HERDR, "pane", "rename", roles[role], PANE_LABELS[role]], runner=runner)

    return PaneSet(tab_id=tab_id, roles=roles)


def close_layout(tab_id: str, runner: Runner = _subprocess_runner) -> None:
    herdr_json([HERDR, "tab", "close", tab_id], runner=runner)


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="council panes",
        description="Create or tear down the Herdr pane layout for the council.",
    )
    sub = parser.add_subparsers(dest="action", required=True)

    create = sub.add_parser("create", help="create the council layout")
    create.add_argument("--label", default=DEFAULT_TAB_LABEL)
    create.add_argument("--cwd", default=".")
    create.add_argument("--state", type=Path, help="write the pane set here for later stages")
    create.add_argument("--focus", action="store_true", help="switch to the tab once built")
    create.add_argument(
        "--replace",
        action="store_true",
        help="close an existing tab with the same label first",
    )

    close = sub.add_parser("close", help="close a council layout")
    close.add_argument("--label", default=DEFAULT_TAB_LABEL)

    args = parser.parse_args(argv)

    if args.action == "close":
        tab_id = find_tab_by_label(args.label)
        if tab_id is None:
            print(f"no tab labelled {args.label!r}")
            return 1
        close_layout(tab_id)
        print(f"closed {tab_id}")
        return 0

    if args.replace:
        existing = find_tab_by_label(args.label)
        if existing:
            close_layout(existing)

    panes = create_layout(label=args.label, cwd=args.cwd, focus=args.focus)
    if args.state:
        panes.save(args.state)
    print(f"tab={panes.tab_id}")
    for role in ROLES:
        print(f"  {role:<13} {panes.roles[role]}  {PANE_LABELS[role]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
