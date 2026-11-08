# ABOUTME: Tests that every council prompt prohibits tools before the task data.
# ABOUTME: Pins the escape-hatch exception, which permits source reads and nothing else.

"""Check tool instructions at each council stage."""

import sys
from datetime import date
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.deliberate import build_chair_prompt, build_review_prompt
from council.fanout import build_stage1_prompt
from council.research import build_distill_prompt


@pytest.mark.parametrize("stage", ["answer", "review", "chair", "distil"])
def test_headless_stages_prohibit_tools_before_the_task_data(stage):
    question = "Check whether Friday 11 September 2026 is correct."
    answers = {"one": "Friday 11 September 2026", "two": "Check the calendar."}
    labels = {"one": "A", "two": "B"}
    prompts = {
        "answer": lambda: build_stage1_prompt(question, brief=None),
        "review": lambda: build_review_prompt(question, answers, labels, brief=None),
        "chair": lambda: build_chair_prompt(
            question, answers, labels, order=["A", "B"], brief=None, rankings={}
        ),
        "distil": lambda: build_distill_prompt(question, "[S1] Friday", date(2026, 9, 7)),
    }
    prompt = prompts[stage]()
    instructions = prompt[:prompt.index(question)]
    assert "Do not call tools." in instructions
    assert "Do not run shell commands" in instructions
    assert "state the uncertainty" in instructions
    if stage == "review":
        assert "Do not carry out the original request yourself." in instructions


def test_escape_hatch_permits_only_source_reads():
    prompt = build_stage1_prompt("question", brief="[S1] fact", escape_hatch_slots=1)
    assert "Use only the file-read tool for the source paths listed below." in prompt
    assert "Do not run shell commands" in prompt
    assert "Do not call tools." not in prompt
    assert "Use source reads only to check claims against the listed sources." in prompt
    assert "/config/.llm-council/current/raw/S1.md" in prompt


def test_escape_hatch_without_sources_does_not_permit_tools():
    prompt = build_stage1_prompt("question", brief=None, escape_hatch_slots=1)
    assert "Do not call tools." in prompt
