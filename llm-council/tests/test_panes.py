# ABOUTME: Tests for the Herdr pane scaffolding used to host the council.
# ABOUTME: Covers the layout plan, socket-API envelope handling and role/pane bookkeeping.

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.panes import (  # noqa: E402
    LAYOUT_PLAN,
    MEMBER_C,
    ORCHESTRATOR,
    PANE_LABELS,
    ROLES,
    HerdrError,
    PaneSet,
    close_layout,
    create_layout,
    find_tab_by_label,
    herdr_json,
    herdr_ok,
)


class FakeHerd:
    """Stands in for the Herdr socket API, recording every command issued."""

    def __init__(self, tabs=None):
        self.commands = []
        self._next_pane = 100
        self._tabs = tabs or []

    def __call__(self, cmd, timeout):
        self.commands.append(list(cmd))
        verb = cmd[1:3]

        if verb == ["tab", "create"]:
            self._next_pane += 1
            return 0, json.dumps(
                {
                    "result": {
                        "tab": {"tab_id": "w9:t1"},
                        "root_pane": {"pane_id": f"w9:p{self._next_pane}"},
                    }
                },
            ), ""

        if verb == ["pane", "split"]:
            self._next_pane += 1
            return 0, json.dumps({"result": {"pane": {"pane_id": f"w9:p{self._next_pane}"}}}), ""

        if verb == ["tab", "list"]:
            return 0, json.dumps({"result": {"tabs": self._tabs}}), ""

        return 0, json.dumps({"result": {"type": "ok"}}), ""

    def issued(self, *prefix):
        return [c for c in self.commands if c[1 : 1 + len(prefix)] == list(prefix)]


# --- layout plan (pure) ---------------------------------------------------


def test_layout_plan_produces_every_role_exactly_once():
    produced = [step.new_role for step in LAYOUT_PLAN]
    assert sorted(produced + [ORCHESTRATOR]) == sorted(ROLES)


def test_layout_plan_only_splits_panes_that_already_exist():
    """A split cannot target a pane the plan has not created yet."""
    existing = {ORCHESTRATOR}
    for step in LAYOUT_PLAN:
        assert step.source_role in existing, f"{step.source_role} split before it exists"
        existing.add(step.new_role)


def test_layout_plan_ratios_are_fractions():
    for step in LAYOUT_PLAN:
        assert 0.0 < step.ratio < 1.0
        assert step.direction in ("right", "down")


def test_every_role_has_a_human_label():
    assert set(PANE_LABELS) == set(ROLES)


def test_member_c_label_names_its_model():
    assert PANE_LABELS[MEMBER_C] == "C · gemini 4 argon"


# --- herdr_json -----------------------------------------------------------


def test_herdr_json_unwraps_the_result_envelope():
    def runner(cmd, timeout):
        return 0, json.dumps({"id": "cli:x", "result": {"type": "ok", "n": 1}}), ""

    assert herdr_json(["herdr", "pane", "list"], runner=runner) == {"type": "ok", "n": 1}


def test_herdr_json_raises_on_nonzero_exit():
    def runner(cmd, timeout):
        return 1, "", "socket missing"

    with pytest.raises(HerdrError, match="socket missing"):
        herdr_json(["herdr", "pane", "list"], runner=runner)


def test_herdr_json_raises_on_unparseable_output():
    def runner(cmd, timeout):
        return 0, "not json at all", ""

    with pytest.raises(HerdrError):
        herdr_json(["herdr", "pane", "list"], runner=runner)


def test_herdr_json_raises_when_envelope_has_no_result():
    def runner(cmd, timeout):
        return 0, json.dumps({"id": "cli:x", "error": "nope"}), ""

    with pytest.raises(HerdrError):
        herdr_json(["herdr", "pane", "list"], runner=runner)


def test_herdr_ok_accepts_empty_output():
    """`herdr pane run` succeeds silently; demanding JSON from it would be a false failure."""

    def runner(cmd, timeout):
        return 0, "", ""

    herdr_ok(["herdr", "pane", "run", "w1:p1", "bash", "x.sh"], runner=runner)


def test_herdr_ok_raises_on_nonzero_exit():
    def runner(cmd, timeout):
        return 1, "", "no such pane"

    with pytest.raises(HerdrError, match="no such pane"):
        herdr_ok(["herdr", "pane", "run", "w1:p1", "bash", "x.sh"], runner=runner)


# --- create_layout --------------------------------------------------------


def test_create_layout_returns_a_pane_for_every_role():
    herd = FakeHerd()
    panes = create_layout(label="council", cwd="/tmp/run", runner=herd)
    assert set(panes.roles) == set(ROLES)
    assert panes.tab_id == "w9:t1"
    assert len(set(panes.roles.values())) == len(ROLES), "roles must map to distinct panes"


def test_create_layout_creates_the_tab_before_splitting():
    herd = FakeHerd()
    create_layout(label="council", cwd="/tmp/run", runner=herd)
    verbs = [c[1:3] for c in herd.commands]
    assert verbs[0] == ["tab", "create"]
    assert ["pane", "split"] in verbs


def test_create_layout_issues_one_split_per_plan_step():
    herd = FakeHerd()
    create_layout(label="council", cwd="/tmp/run", runner=herd)
    assert len(herd.issued("pane", "split")) == len(LAYOUT_PLAN)


def test_create_layout_passes_cwd_so_panes_start_in_the_run_directory():
    herd = FakeHerd()
    create_layout(label="council", cwd="/tmp/run", runner=herd)
    create_cmd = herd.issued("tab", "create")[0]
    assert "--cwd" in create_cmd
    assert "/tmp/run" in create_cmd
    for split in herd.issued("pane", "split"):
        assert "--cwd" in split


def test_create_layout_does_not_steal_focus():
    """The council runs in the background; it must not yank the user out of their tab."""
    herd = FakeHerd()
    create_layout(label="council", cwd="/tmp/run", runner=herd, focus=False)
    for cmd in herd.issued("tab", "create") + herd.issued("pane", "split"):
        assert "--no-focus" in cmd


def test_create_layout_labels_every_pane():
    herd = FakeHerd()
    create_layout(label="council", cwd="/tmp/run", runner=herd)
    renames = herd.issued("pane", "rename")
    assert len(renames) == len(ROLES)


def test_create_layout_splits_the_pane_named_by_the_plan():
    herd = FakeHerd()
    panes = create_layout(label="council", cwd="/tmp/run", runner=herd)
    splits = herd.issued("pane", "split")
    for step, cmd in zip(LAYOUT_PLAN, splits):
        assert panes.roles[step.source_role] in cmd
        assert "--direction" in cmd and step.direction in cmd


# --- find_tab_by_label ----------------------------------------------------


def test_find_tab_by_label_returns_matching_tab_id():
    herd = FakeHerd(tabs=[{"tab_id": "w1:t3", "label": "council"}])
    assert find_tab_by_label("council", runner=herd) == "w1:t3"


def test_find_tab_by_label_returns_none_when_absent():
    herd = FakeHerd(tabs=[{"tab_id": "w1:t1", "label": "Claude"}])
    assert find_tab_by_label("council", runner=herd) is None


# --- close_layout ---------------------------------------------------------


def test_close_layout_closes_the_tab():
    herd = FakeHerd()
    close_layout("w9:t1", runner=herd)
    assert herd.issued("tab", "close")[0][-1] == "w9:t1"


# --- PaneSet persistence --------------------------------------------------


def test_paneset_round_trips_through_disk(tmp_path):
    """Later stages re-read the layout rather than rebuilding it."""
    original = PaneSet(tab_id="w9:t1", roles={role: f"w9:p{i}" for i, role in enumerate(ROLES)})
    path = tmp_path / "panes.json"
    original.save(path)
    assert PaneSet.load(path) == original


def test_paneset_rejects_a_missing_role():
    with pytest.raises(ValueError):
        PaneSet(tab_id="w9:t1", roles={ORCHESTRATOR: "w9:p1"})
