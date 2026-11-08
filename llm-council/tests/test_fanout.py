# ABOUTME: Tests for stage 1 fan-out: prompt assembly, per-member scripts and sentinel waiting.
# ABOUTME: Verifies dispatch never depends on scraping terminal output.

import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.fanout import (  # noqa: E402
    AGY_EFFORT,
    CHAIRMAN_MEMBER,
    DEBATE_STANCES,
    HATCH_CAPABLE_KINDS,
    MEMBERS,
    AgentFailedError,
    Member,
    PromptTooLargeError,
    StageTimeout,
    build_stage1_prompt,
    collect_answers,
    detect_hatch_usage,
    dispatch_stage,
    member_paths,
    render_member_script,
    wait_for_stage,
)
from council.panes import MEMBER_A, MEMBER_B, MEMBER_C, ROLES  # noqa: E402


def paths_for(tmp_path, member, stage="stage1"):
    return member_paths(tmp_path, member.role, stage)


# --- roster ---------------------------------------------------------------


def test_roster_has_three_distinct_members_from_three_vendors():
    assert len(MEMBERS) == 3
    assert len({m.role for m in MEMBERS}) == 3
    assert {m.kind for m in MEMBERS} == {"claude", "codex", "agy"}


def test_roster_roles_are_real_pane_roles():
    for member in MEMBERS:
        assert member.role in ROLES


def test_roster_labels_are_anonymous_single_letters():
    """Members are shown to each other as A/B/C, never by vendor."""
    assert sorted(m.label for m in MEMBERS) == ["A", "B", "C"]


# --- paths ----------------------------------------------------------------


def test_member_paths_are_scoped_to_stage_and_role(tmp_path):
    p = member_paths(tmp_path, MEMBER_A, "stage1")
    assert p.answer.parent == tmp_path / "stage1"
    assert MEMBER_A in p.answer.name
    assert p.done.suffix == ".done"
    assert p.status.suffix == ".status"


# --- prompt ---------------------------------------------------------------


def test_stage1_prompt_contains_the_question():
    prompt = build_stage1_prompt("why is the sky blue?", brief=None)
    assert "why is the sky blue?" in prompt


def test_stage1_prompt_embeds_the_brief_when_present():
    prompt = build_stage1_prompt("q", brief="# Evidence brief\n[S1] fact")
    assert "[S1] fact" in prompt


def test_stage1_prompt_without_a_brief_does_not_reference_evidence():
    prompt = build_stage1_prompt("q", brief=None)
    assert "evidence brief" not in prompt.lower()


def test_stage1_prompt_asks_for_a_bare_answer():
    """Answers are compared side by side, so preamble is noise."""
    prompt = build_stage1_prompt("q", brief=None)
    assert "preamble" in prompt.lower()


def test_brief_wins_over_recollection_on_a_direct_contradiction():
    """Retrieval is recency-filtered; training data is the older source."""
    prompt = build_stage1_prompt("q", brief="[S1] fact").lower()
    assert "contradict" in prompt or "disagree" in prompt
    assert "brief" in prompt


def test_prior_knowledge_is_permitted_where_the_brief_is_thin():
    """Five scraped sources are not a survey; suppressing recall would lose real coverage."""
    prompt = build_stage1_prompt("q", brief="[S1] fact").lower()
    assert "silent" in prompt or "thin" in prompt or "does not cover" in prompt


def test_prior_knowledge_claims_must_be_marked():
    prompt = build_stage1_prompt("q", brief="[S1] fact")
    assert "[prior]" in prompt


def test_grounding_rules_are_absent_without_a_brief():
    """With no brief there is nothing to defer to, so the rules would be noise."""
    prompt = build_stage1_prompt("q", brief=None)
    assert "[prior]" not in prompt


def test_escape_hatch_lists_literal_paths_a_member_can_open():
    prompt = build_stage1_prompt("q", brief="[S1] fact", escape_hatch_slots=2)
    assert "/config/.llm-council/current/raw/S1.md" in prompt
    assert "/config/.llm-council/current/raw/S2.md" in prompt


def test_escape_hatch_is_absent_by_default():
    """Members must have symmetric access; the hatch is opt-in once all three can read."""
    prompt = build_stage1_prompt("q", brief="[S1] fact")
    assert "current/raw" not in prompt


def test_escape_hatch_needs_a_brief_to_be_meaningful():
    prompt = build_stage1_prompt("q", brief=None, escape_hatch_slots=3)
    assert "current/raw" not in prompt


def test_codex_gets_no_tools_even_with_the_hatch_on(tmp_path):
    """Codex's sandbox fails closed here, so granting it the hatch would be a lie."""
    member = next(m for m in MEMBERS if m.kind == "codex")
    script = render_member_script(
        member, paths_for(tmp_path, member), timeout_s=900, escape_hatch=True
    )
    assert "--sandbox read-only" in script
    assert "danger-full-access" not in script


def test_hatch_capable_kinds_excludes_codex():
    assert "codex" not in HATCH_CAPABLE_KINDS
    assert {"claude", "agy"} <= set(HATCH_CAPABLE_KINDS)


def test_hatch_prompt_asks_for_a_usage_marker():
    prompt = build_stage1_prompt("q", brief="[S1] fact", escape_hatch_slots=2)
    assert "[Sn:full]" in prompt


def test_detect_hatch_usage_reports_sources_opened_per_member():
    usage = detect_hatch_usage(
        {
            MEMBER_A: "Per [S1:full] the flag was renamed, though [S2] disagrees.",
            MEMBER_B: "Only the brief here [S2].",
            MEMBER_C: "[S3:full] and again [S3:full] and [S1:full].",
        }
    )
    assert usage[MEMBER_A] == [1]
    assert usage[MEMBER_B] == []
    assert usage[MEMBER_C] == [1, 3], "duplicates collapse, order is stable"


def test_debate_stances_cover_every_member_and_differ():
    assert set(DEBATE_STANCES) == {m.role for m in MEMBERS}
    assert len(set(DEBATE_STANCES.values())) == 3


def test_stage1_prompt_embeds_a_stance_when_given():
    prompt = build_stage1_prompt("q", brief=None, stance="Argue the strongest case.")
    assert "Argue the strongest case." in prompt


def test_stage1_prompt_omits_stance_language_by_default():
    """Default mode compares models; stance would confound that."""
    default = build_stage1_prompt("q", brief=None)
    stanced = build_stage1_prompt("q", brief=None, stance="Argue the strongest case.")
    assert len(default) < len(stanced)


# --- script rendering -----------------------------------------------------


def test_claude_script_pipes_prompt_via_stdin_and_disables_tools(tmp_path):
    member = next(m for m in MEMBERS if m.kind == "claude")
    script = render_member_script(member, paths_for(tmp_path, member), timeout_s=900)
    assert "claude -p" in script
    assert "--model claude-opus-5-5" in script
    assert '--allowedTools ""' in script
    assert "< " in script, "prompt should arrive on stdin, not argv"


def test_claude_gets_read_only_when_the_escape_hatch_is_on(tmp_path):
    member = next(m for m in MEMBERS if m.kind == "claude")
    script = render_member_script(
        member, paths_for(tmp_path, member), timeout_s=900, escape_hatch=True
    )
    assert '--allowedTools "Read"' in script
    # Read only: a member must never edit or run commands.
    for forbidden in ("Edit", "Write", "Bash"):
        assert forbidden not in script


def test_chairman_script_pins_fable_because_synthesis_gets_the_strongest_model(tmp_path):
    """The chair reads every answer and every review; it gets the most capable model."""
    script = render_member_script(
        CHAIRMAN_MEMBER, paths_for(tmp_path, CHAIRMAN_MEMBER, "chair"), timeout_s=900
    )
    assert "claude -p" in script
    assert "--model claude-fable-5-1" in script


def test_codex_script_captures_via_output_file_not_stdout(tmp_path):
    """codex prints token accounting to stdout, so -o is the only clean capture."""
    member = next(m for m in MEMBERS if m.kind == "codex")
    p = paths_for(tmp_path, member)
    script = render_member_script(member, p, timeout_s=900)
    assert "codex exec -" in script
    assert "--skip-git-repo-check" in script
    assert "--sandbox read-only" in script
    assert "--model gpt-6-astra" in script
    assert f"-o {p.answer}" in script


def test_agy_script_uses_argv_because_stdin_is_ignored(tmp_path):
    member = next(m for m in MEMBERS if m.kind == "agy")
    p = paths_for(tmp_path, member)
    script = render_member_script(member, p, timeout_s=900)
    assert "agy -p" in script
    assert "$(cat" in script, "agy ignores stdin; prompt must be substituted into argv"
    assert "--model gemini-4-argon" in script
    assert "--print-timeout" in script


@pytest.mark.parametrize("member", MEMBERS, ids=lambda m: m.kind)
def test_every_script_streams_to_the_pane_as_well_as_to_disk(tmp_path, member):
    """The pane is the live view: fully redirecting output would leave it blank."""
    p = paths_for(tmp_path, member)
    script = render_member_script(member, p, timeout_s=900)
    assert "| tee" in script, "output must reach the terminal, not only a file"
    # codex writes its answer with -o, so its visible stream is the progress log.
    tee_target = p.log if member.kind == "codex" else p.answer
    assert f"tee {tee_target}" in script or f"tee '{tee_target}'" in script


@pytest.mark.parametrize("member", MEMBERS, ids=lambda m: m.kind)
def test_every_script_records_status_then_signals_completion(tmp_path, member):
    p = paths_for(tmp_path, member)
    script = render_member_script(member, p, timeout_s=900)
    assert str(p.status) in script
    assert str(p.done) in script
    assert script.index(str(p.status)) < script.rindex(str(p.done))


@pytest.mark.parametrize("member", MEMBERS, ids=lambda m: m.kind)
def test_no_script_uses_errexit(tmp_path, member):
    """`set -e` would abort before the sentinel, hanging the orchestrator forever."""
    script = render_member_script(member, paths_for(tmp_path, member), timeout_s=900)
    directives = [
        line.strip()
        for line in script.splitlines()
        if line.strip().startswith("set ")  # ignore comments discussing errexit
    ]
    for directive in directives:
        flags = directive.removeprefix("set ").split()
        assert not any(f.startswith("-") and "e" in f.lstrip("-") for f in flags), directive


def test_script_signals_completion_even_when_the_agent_fails(tmp_path, monkeypatch):
    """Executed for real against a stub binary that exits non-zero."""
    member = next(m for m in MEMBERS if m.kind == "claude")
    p = paths_for(tmp_path, member)
    p.prompt.parent.mkdir(parents=True, exist_ok=True)
    p.prompt.write_text("question")

    stub_dir = tmp_path / "stub"
    stub_dir.mkdir()
    stub = stub_dir / "claude"
    stub.write_text("#!/usr/bin/env bash\necho 'boom' >&2\nexit 3\n")
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)

    script_path = tmp_path / "run.sh"
    script_path.write_text(render_member_script(member, p, timeout_s=30))
    env = dict(os.environ, PATH=f"{stub_dir}:{os.environ['PATH']}")
    subprocess.run(["bash", str(script_path)], env=env, timeout=60, check=False)

    assert p.done.exists(), "sentinel must fire even on failure"
    assert p.status.read_text().strip() == "3"
    assert "boom" in p.log.read_text()


def test_agy_prompt_over_the_argv_limit_is_rejected(tmp_path):
    member = next(m for m in MEMBERS if m.kind == "agy")
    p = paths_for(tmp_path, member)
    p.prompt.parent.mkdir(parents=True, exist_ok=True)
    p.prompt.write_text("z" * 200_000)
    with pytest.raises(PromptTooLargeError):
        render_member_script(member, p, timeout_s=900)


def test_large_prompt_is_fine_for_stdin_agents(tmp_path):
    """Only agy is argv-bound; claude and codex read stdin and have no such cap."""
    member = next(m for m in MEMBERS if m.kind == "claude")
    p = paths_for(tmp_path, member)
    p.prompt.parent.mkdir(parents=True, exist_ok=True)
    p.prompt.write_text("z" * 200_000)
    render_member_script(member, p, timeout_s=900)


# --- dispatch -------------------------------------------------------------


class FakeHerd:
    def __init__(self):
        self.commands = []

    def __call__(self, cmd, timeout):
        self.commands.append(list(cmd))
        return 0, json.dumps({"result": {"type": "ok"}}), ""


def test_dispatch_runs_each_member_in_its_own_pane(tmp_path):
    herd = FakeHerd()
    panes = {MEMBER_A: "w1:p1", MEMBER_B: "w1:p2", MEMBER_C: "w1:p3"}
    dispatch_stage(
        panes=panes,
        run_dir=tmp_path,
        stage="stage1",
        prompts={m.role: "the prompt" for m in MEMBERS},
        runner=herd,
    )
    assert len(herd.commands) == 3
    dispatched_panes = {c[3] for c in herd.commands}
    assert dispatched_panes == set(panes.values())
    for cmd in herd.commands:
        assert cmd[1:3] == ["pane", "run"]
        assert "bash" in cmd


def test_dispatch_writes_prompt_and_script_to_disk(tmp_path):
    dispatch_stage(
        panes={m.role: f"w1:p{i}" for i, m in enumerate(MEMBERS)},
        run_dir=tmp_path,
        stage="stage1",
        prompts={m.role: f"prompt for {m.label}" for m in MEMBERS},
        runner=FakeHerd(),
    )
    for member in MEMBERS:
        p = member_paths(tmp_path, member.role, "stage1")
        assert p.prompt.read_text() == f"prompt for {member.label}"
        assert p.script.exists()


def test_dispatch_clears_stale_sentinels_from_a_previous_run(tmp_path):
    """A leftover .done would make the orchestrator read the previous answer."""
    stale = member_paths(tmp_path, MEMBER_A, "stage1")
    stale.done.parent.mkdir(parents=True, exist_ok=True)
    stale.done.write_text("")
    stale.status.write_text("0")

    dispatch_stage(
        panes={m.role: "w1:p1" for m in MEMBERS},
        run_dir=tmp_path,
        stage="stage1",
        prompts={m.role: "p" for m in MEMBERS},
        runner=FakeHerd(),
    )
    assert not stale.done.exists()
    assert not stale.status.exists()


# --- waiting --------------------------------------------------------------


def test_wait_returns_exit_status_once_every_sentinel_appears(tmp_path):
    roles = [m.role for m in MEMBERS]
    for role in roles:
        p = member_paths(tmp_path, role, "stage1")
        p.done.parent.mkdir(parents=True, exist_ok=True)
        p.status.write_text("0")
        p.done.write_text("")

    result = wait_for_stage(roles, tmp_path, "stage1", timeout=5, poll=0.01)
    assert result == {role: 0 for role in roles}


def test_wait_times_out_when_a_member_never_finishes(tmp_path):
    roles = [m.role for m in MEMBERS]
    done = member_paths(tmp_path, roles[0], "stage1")
    done.done.parent.mkdir(parents=True, exist_ok=True)
    done.status.write_text("0")
    done.done.write_text("")

    with pytest.raises(StageTimeout) as excinfo:
        wait_for_stage(roles, tmp_path, "stage1", timeout=0.05, poll=0.01)
    # The error must name who stalled, not just that something did.
    assert roles[1] in str(excinfo.value)


# --- collection -----------------------------------------------------------


def test_collect_returns_each_members_answer(tmp_path):
    roles = [m.role for m in MEMBERS]
    for i, role in enumerate(roles):
        p = member_paths(tmp_path, role, "stage1")
        p.answer.parent.mkdir(parents=True, exist_ok=True)
        p.answer.write_text(f"answer {i}")
        p.status.write_text("0")
    answers = collect_answers(roles, tmp_path, "stage1")
    assert answers[roles[0]] == "answer 0"


def test_collect_raises_when_a_member_exited_non_zero(tmp_path):
    roles = [m.role for m in MEMBERS]
    for role in roles:
        p = member_paths(tmp_path, role, "stage1")
        p.answer.parent.mkdir(parents=True, exist_ok=True)
        p.answer.write_text("ok")
        p.status.write_text("0")
    bad = member_paths(tmp_path, roles[2], "stage1")
    bad.status.write_text("1")
    bad.log.write_text("rate limited")

    with pytest.raises(AgentFailedError, match="rate limited"):
        collect_answers(roles, tmp_path, "stage1")


def test_collect_raises_on_empty_answer(tmp_path):
    """agy exits 0 with no output when a permission is auto-denied."""
    roles = [m.role for m in MEMBERS]
    for role in roles:
        p = member_paths(tmp_path, role, "stage1")
        p.answer.parent.mkdir(parents=True, exist_ok=True)
        p.answer.write_text("ok")
        p.status.write_text("0")
    member_paths(tmp_path, roles[1], "stage1").answer.write_text("  \n ")

    with pytest.raises(AgentFailedError):
        collect_answers(roles, tmp_path, "stage1")


def test_agy_script_passes_effort_because_the_model_requires_it(tmp_path):
    """Member C's effort is pinned, not left to agy's default."""
    member = next(m for m in MEMBERS if m.kind == "agy")
    p = paths_for(tmp_path, member)
    script = render_member_script(member, p, timeout_s=900)
    assert f"--effort {AGY_EFFORT}" in script
    assert AGY_EFFORT in {"low", "medium", "high"}
