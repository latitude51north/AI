# ABOUTME: Tests for end-to-end council orchestration, helpers and the full pipeline.
# ABOUTME: Covers run-directory naming, ballot tolerance, the final document and stage sequencing.

import json
import sys
from datetime import datetime
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.panes import CHAIRMAN, MEMBER_A, MEMBER_B, MEMBER_C  # noqa: E402
from council.run import (  # noqa: E402
    build_final_document,
    default_run_dir,
    main,
    parse_ballots,
    run_council,
)

ROLES = [MEMBER_A, MEMBER_B, MEMBER_C]
LABELS = {MEMBER_A: "C", MEMBER_B: "A", MEMBER_C: "B"}


# --- run directory --------------------------------------------------------


def test_default_run_dir_is_timestamped_and_readable():
    path = default_run_dir("Why is the sky blue?", now=datetime(2026, 7, 27, 14, 30, 5))
    assert "20260727-143005" in path.name
    assert "why-is-the-sky-blue" in path.name


def test_default_run_dir_is_bounded_for_a_long_question():
    path = default_run_dir("word " * 200, now=datetime(2026, 7, 27, 14, 30, 5))
    assert len(path.name) < 100


def test_default_run_dir_survives_a_question_with_no_word_characters():
    path = default_run_dir("???", now=datetime(2026, 7, 27, 14, 30, 5))
    assert path.name.endswith("query")


# --- argument guards ------------------------------------------------------


def test_query_without_web_is_rejected_rather_than_ignored(capsys):
    """Retrieval never runs without --web, so a search query would silently vanish."""
    with pytest.raises(SystemExit) as exc:
        main(["a question", "--query", "short searchable terms"])
    assert exc.value.code == 2
    assert "--query only applies with --web" in capsys.readouterr().err


# --- ballots --------------------------------------------------------------


def test_parse_ballots_keeps_the_good_and_reports_the_bad():
    reviews = {
        MEMBER_A: "reasons\nRANKING: A > B > C",
        MEMBER_B: "I liked them all equally.",
        MEMBER_C: "RANKING: B > A > C",
    }
    ballots, failures = parse_ballots(reviews, expected={"A", "B", "C"}, labels=LABELS)
    assert set(ballots) == {"C", "B"}  # keyed by reviewer label, not role
    assert ballots["C"] == ["A", "B", "C"]
    assert MEMBER_B in failures


def test_parse_ballots_returns_empty_when_every_reviewer_fails():
    """A total ballot failure must degrade, not crash: the chair can still synthesise."""
    reviews = {role: "no ranking here" for role in ROLES}
    ballots, failures = parse_ballots(reviews, expected={"A", "B", "C"}, labels=LABELS)
    assert ballots == {}
    assert len(failures) == 3


# --- final document -------------------------------------------------------


def test_final_document_leads_with_the_chair_answer():
    doc = build_final_document(
        question="why?",
        chair_answer="Because of scattering.",
        order=["A", "B", "C"],
        divergence={"mean_similarity": 0.3, "converged": False, "pairs": []},
        labels=LABELS,
        members={MEMBER_A: "claude opus", MEMBER_B: "codex", MEMBER_C: "gemini"},
    )
    assert doc.index("Because of scattering.") < doc.index("mean similarity")


def test_final_document_reveals_the_mapping_to_the_reader():
    """Anonymity protects the reviewers, not the user: they should see who was who."""
    doc = build_final_document(
        question="why?",
        chair_answer="answer",
        order=["A", "B", "C"],
        divergence={"mean_similarity": 0.3, "converged": False, "pairs": []},
        labels=LABELS,
        members={MEMBER_A: "claude opus", MEMBER_B: "codex", MEMBER_C: "gemini"},
    )
    assert "claude opus" in doc
    assert "codex" in doc
    # member_a was shown to the panel as C
    assert "C" in doc


def test_final_document_warns_when_members_converged():
    doc = build_final_document(
        question="why?",
        chair_answer="answer",
        order=["A", "B", "C"],
        divergence={"mean_similarity": 0.95, "converged": True, "pairs": []},
        labels=LABELS,
        members={r: "m" for r in ROLES},
    )
    assert "converged" in doc.lower()


def test_final_document_handles_a_missing_panel_ranking():
    doc = build_final_document(
        question="why?",
        chair_answer="answer",
        order=[],
        divergence={"mean_similarity": 0.3, "converged": False, "pairs": []},
        labels=LABELS,
        members={r: "m" for r in ROLES},
    )
    assert "answer" in doc


# --- full pipeline (integration, every subprocess boundary faked) ---------


class FakeCouncil:
    """Stands in for every process the council shells out to.

    Routes on the executable name, so one fake covers Herdr, Firecrawl, agy and
    uv. Crucially it also plays the part of the *pane*: when Herdr is told to run
    a member's script, this writes the answer, status and sentinel that the real
    script would have written, so stage sequencing is exercised for real rather
    than stubbed out.
    """

    def __init__(self, answers=None, reviews=None, chair="The synthesised answer.", fail=()):
        self.answers = answers or {}
        self.reviews = reviews or {}
        self.chair = chair
        self.fail = set(fail)
        self.commands: list[list[str]] = []
        self.stages: list[str] = []
        self.panes_run: list[str] = []
        self._next_pane = 0

    # -- what a member's pane would have produced --------------------------
    def _output_for(self, stage: str, role: str) -> str:
        if stage == "stage1":
            return self.answers.get(role, f"Distinct answer from {role} about topic {role}.")
        if stage == "stage2":
            return self.reviews.get(role, "Justification here.\n\nRANKING: A > B > C")
        return self.chair

    def _run_script(self, script: str) -> None:
        path = Path(script)
        role, stage = path.stem, path.parent.name
        if stage not in self.stages:
            self.stages.append(stage)
        self.panes_run.append(f"{stage}/{role}")

        base = path.parent
        if role in self.fail:
            (base / f"{role}.status").write_text("1\n")
            (base / f"{role}.log").write_text("the agent exploded\n")
        else:
            (base / f"{role}.md").write_text(self._output_for(stage, role))
            (base / f"{role}.status").write_text("0\n")
        (base / f"{role}.done").touch()

    # -- the runner interface ---------------------------------------------
    def __call__(self, cmd, timeout):
        cmd = list(cmd)
        self.commands.append(cmd)
        exe = cmd[0]

        if exe == "herdr":
            return self._herdr(cmd)
        if exe == "firecrawl":
            out = Path(cmd[cmd.index("-o") + 1])
            out.write_text(
                json.dumps(
                    {
                        "success": True,
                        "creditsUsed": 3,
                        "data": {
                            "web": [
                                {"url": "https://a.test", "title": "A", "markdown": "evidence body"}
                            ]
                        },
                    }
                )
            )
            return 0, "", ""
        if exe == "agy":
            return 0, "distilled findings [S1]", ""
        if exe == "uv":
            Path(cmd[-1]).write_text("converted")
            return 0, "", ""
        raise AssertionError(f"unexpected command: {exe}")

    def _herdr(self, cmd):
        action = " ".join(cmd[1:3])
        if action == "tab create":
            return 0, json.dumps({"result": {"tab": {"tab_id": "T1"}, "root_pane": {"pane_id": "p0"}}}), ""
        if action == "pane split":
            self._next_pane += 1
            return 0, json.dumps({"result": {"pane": {"pane_id": f"p{self._next_pane}"}}}), ""
        if action == "pane rename":
            return 0, json.dumps({"result": {}}), ""
        if action == "tab list":
            return 0, json.dumps({"result": {"tabs": [{"label": "council", "tab_id": "T1"}]}}), ""
        if action == "tab close":
            return 0, json.dumps({"result": {}}), ""
        if action == "pane run":
            if cmd[4] == "bash":
                self._run_script(cmd[5])
            return 0, "", ""
        raise AssertionError(f"unexpected herdr command: {action}")


def test_run_council_completes_and_writes_the_final_document(tmp_path):
    fake = FakeCouncil()
    final = run_council("why is the sky blue?", run_dir=tmp_path, runner=fake)

    assert final == tmp_path / "final.md"
    text = final.read_text()
    assert "The synthesised answer." in text
    assert "why is the sky blue?" in text


def test_run_council_runs_the_three_stages_in_order(tmp_path):
    fake = FakeCouncil()
    run_council("q", run_dir=tmp_path, runner=fake)
    assert fake.stages == ["stage1", "stage2", "stage3"]


def test_run_council_asks_every_member_then_only_the_chair(tmp_path):
    """Stages 1 and 2 use all three members; stage 3 is the chair alone."""
    fake = FakeCouncil()
    run_council("q", run_dir=tmp_path, runner=fake)

    assert sorted(p for p in fake.panes_run if p.startswith("stage1")) == [
        f"stage1/{r}" for r in sorted(ROLES)
    ]
    assert [p for p in fake.panes_run if p.startswith("stage3")] == [f"stage3/{CHAIRMAN}"]


def test_run_council_records_the_full_audit_trail(tmp_path):
    fake = FakeCouncil()
    run_council("q", run_dir=tmp_path, runner=fake)

    for name in ("query.txt", "labels.json", "rankings.json", "divergence.json", "panes.json"):
        assert (tmp_path / name).exists(), f"{name} was not written"


def test_run_council_survives_every_ballot_being_unusable(tmp_path):
    """One bad stage must not sink a run that otherwise produced answers."""
    fake = FakeCouncil(reviews={r: "I decline to rank these." for r in ROLES})
    final = run_council("q", run_dir=tmp_path, runner=fake)

    rankings = json.loads((tmp_path / "rankings.json").read_text())
    assert rankings["ballots"] == {}
    assert len(rankings["failures"]) == 3
    assert rankings["order"] == ["A", "B", "C"]  # fallback, not a real verdict
    assert final.exists()


def test_run_council_reports_a_member_that_failed(tmp_path):
    fake = FakeCouncil(fail={MEMBER_B})
    with pytest.raises(Exception) as exc:
        run_council("q", run_dir=tmp_path, runner=fake)
    assert MEMBER_B in str(exc.value)
    assert "exploded" in str(exc.value)


def test_run_council_times_out_naming_the_member_that_stalled(tmp_path):
    class SilentPanes(FakeCouncil):
        def _run_script(self, script):
            if Path(script).stem == MEMBER_C:
                return  # this member never signals
            super()._run_script(script)

    fake = SilentPanes()
    with pytest.raises(Exception) as exc:
        run_council("q", run_dir=tmp_path, runner=fake, timeout=0)
    assert MEMBER_C in str(exc.value)


def test_run_council_retrieves_a_brief_and_shows_it_to_members(tmp_path):
    fake = FakeCouncil()
    run_council("q", run_dir=tmp_path, web=True, limit=1, runner=fake, mirror_dir=tmp_path / "mir")

    assert (tmp_path / "research" / "brief.md").exists()
    prompt = (tmp_path / "stage1" / f"{MEMBER_A}.prompt.md").read_text()
    assert "EVIDENCE BRIEF" in prompt
    assert "distilled findings [S1]" in prompt


def test_run_council_searches_on_the_search_query_when_given(tmp_path):
    fake = FakeCouncil()
    run_council(
        "a very long prompt " * 100,
        run_dir=tmp_path,
        web=True,
        limit=1,
        search_query="short terms",
        runner=fake,
        mirror_dir=tmp_path / "mir",
    )
    searches = [c for c in fake.commands if c[0] == "firecrawl"]
    assert searches[0][2] == "short terms"


def test_run_council_exports_the_synthesis(tmp_path):
    fake = FakeCouncil()
    target = tmp_path / "out" / "synthesis.md"
    run_council("q", run_dir=tmp_path, export=target, runner=fake)
    assert target.read_text() == (tmp_path / "final.md").read_text()


def test_failed_export_does_not_discard_a_completed_run(tmp_path, capsys):
    """Throwing away a finished council over a format conversion would be absurd."""

    class BadExport(FakeCouncil):
        def __call__(self, cmd, timeout):
            if cmd[0] == "uv":
                return 1, "", "pandoc fell over"
            return super().__call__(cmd, timeout)

    fake = BadExport()
    final = run_council("q", run_dir=tmp_path, export=tmp_path / "out.docx", runner=fake)

    assert final.exists()
    assert "The synthesised answer." in final.read_text()
    assert "export failed" in capsys.readouterr().err


def test_run_council_writes_the_synthesis_into_the_vault(tmp_path, isolated_vault, capsys):
    fake = FakeCouncil()
    final = run_council("why is the sky blue?", run_dir=tmp_path, runner=fake)

    notes = list(isolated_vault.glob("*.md"))
    assert len(notes) == 1
    note = notes[0].read_text()
    assert notes[0].name.endswith(" why is the sky blue.md")
    assert note.startswith("---\n")
    assert 'chair: "claude claude-fable-5-1"' in note
    assert final.read_text() in note
    assert f"- Run directory: `{tmp_path}`" in note
    assert f"note: {notes[0]}" in capsys.readouterr().out


def test_vault_note_ranks_members_by_name_not_panel_letter(tmp_path, isolated_vault):
    fake = FakeCouncil()
    run_council("q", run_dir=tmp_path, seed=7, runner=fake)

    rankings = json.loads((tmp_path / "rankings.json").read_text())
    labels = json.loads((tmp_path / "labels.json").read_text())
    role_by_label = {label: role for role, label in labels.items()}
    names = {MEMBER_A: "claude claude-opus-5-5", MEMBER_B: "codex gpt-6-astra",
             MEMBER_C: "agy gemini-4-argon"}
    expected = "ranking:\n" + "".join(
        f'  - "{names[role_by_label[label]]}"\n' for label in rankings["order"]
    )
    note = next(isolated_vault.glob("*.md")).read_text()
    assert expected in note


def test_vault_note_records_no_ranking_when_every_ballot_failed(tmp_path, isolated_vault):
    # The alphabetical fallback order is not a verdict, so it must not look like one.
    fake = FakeCouncil(reviews={r: "I decline to rank these." for r in ROLES})
    run_council("q", run_dir=tmp_path, runner=fake)
    note = next(isolated_vault.glob("*.md")).read_text()
    assert "ranking: []\n" in note


def test_run_council_writes_to_an_explicit_vault_dir(tmp_path):
    fake = FakeCouncil()
    run_council("q", run_dir=tmp_path, vault_dir=tmp_path / "elsewhere", runner=fake)
    assert len(list((tmp_path / "elsewhere").glob("*.md"))) == 1


def test_missing_vault_does_not_discard_a_completed_run(tmp_path, capsys):
    fake = FakeCouncil()
    final = run_council(
        "q", run_dir=tmp_path, vault_dir=tmp_path / "unmounted" / "LLM Council", runner=fake
    )

    assert final.exists()
    err = capsys.readouterr().err
    assert "vault note failed" in err
    assert str(final) in err


def test_run_council_records_hatch_usage_when_the_hatch_is_open(tmp_path):
    fake = FakeCouncil(
        answers={
            MEMBER_A: "I read the source [S1:full] and it says so.",
            MEMBER_B: "Answer from the brief alone.",
            MEMBER_C: "Consulted [S2:full] as well.",
        }
    )
    run_council(
        "q", run_dir=tmp_path, web=True, limit=1, escape_hatch=True,
        runner=fake, mirror_dir=tmp_path / "mir",
    )
    usage = json.loads((tmp_path / "hatch.json").read_text())
    assert usage[MEMBER_A] == [1]
    assert usage[MEMBER_B] == []
    assert usage[MEMBER_C] == [2]


def test_run_council_leaves_the_panes_up_by_default(tmp_path):
    fake = FakeCouncil()
    run_council("q", run_dir=tmp_path, runner=fake)
    assert not any(c[1:3] == ["tab", "close"] for c in fake.commands)


def test_run_council_closes_the_layout_when_asked(tmp_path):
    fake = FakeCouncil()
    run_council("q", run_dir=tmp_path, keep_panes=False, runner=fake)
    assert any(c[1:3] == ["tab", "close"] for c in fake.commands)


def test_run_council_inlines_attachments_into_every_member_prompt(tmp_path):
    """Access must stay symmetric: no member needs file tools to see an attachment."""
    doc = tmp_path / "spec.md"
    doc.write_text("The spec says widgets are blue.")

    fake = FakeCouncil()
    run_council("q", run_dir=tmp_path, attach=[doc], runner=fake)

    for role in ROLES:
        prompt = (tmp_path / "stage1" / f"{role}.prompt.md").read_text()
        assert "[A1] spec.md" in prompt
        assert "widgets are blue" in prompt


def test_run_council_warns_when_the_evidence_is_thinner_than_requested(tmp_path, capsys):
    fake = FakeCouncil()  # its fake search only ever yields one source
    run_council("q", run_dir=tmp_path, web=True, limit=5, runner=fake, mirror_dir=tmp_path / "m")

    out = capsys.readouterr().out
    assert "only 1 of 5 sources" in out


def test_main_reports_a_failure_without_a_traceback(capsys):
    """A bad attachment must surface as a message and exit 1, not a stack trace."""
    code = main(["a question", "--attach", "/nonexistent/nope.md"])
    assert code == 1
    err = capsys.readouterr().err
    assert "council run failed" in err
    assert "nope.md" in err
