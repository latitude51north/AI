# ABOUTME: Pins the argv budget shared by the distiller and member prompts to the real kernel limit.
# ABOUTME: agy takes its prompt through argv, so the ceiling is MAX_ARG_STRLEN rather than a policy choice.

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.fanout import AGY_ARGV_LIMIT  # noqa: E402
from council.research import (  # noqa: E402
    AGY_ARGV_BUDGET,
    DEFAULT_MAX_PROMPT_BYTES,
    MAX_ARG_STRLEN,
)

linux_only = pytest.mark.skipif(
    not sys.platform.startswith("linux"), reason="MAX_ARG_STRLEN is a Linux kernel constant"
)


def _execs_with_argument_of(size: int) -> bool:
    """Can a real exec carry a single argument of this many bytes?"""
    try:
        subprocess.run(["/bin/true", "x" * size], check=True)
    except OSError:  # E2BIG
        return False
    return True


# --- the constant must match reality --------------------------------------


@linux_only
def test_max_arg_strlen_is_the_real_limit_on_this_kernel():
    """The budget is only safe if the constant it derives from is true."""
    assert _execs_with_argument_of(MAX_ARG_STRLEN - 1)
    assert not _execs_with_argument_of(MAX_ARG_STRLEN)


@linux_only
def test_a_prompt_at_the_budget_actually_execs():
    """The whole point of the margin: the budget must be comfortably reachable."""
    assert _execs_with_argument_of(AGY_ARGV_BUDGET)


@linux_only
def test_other_arguments_do_not_count_against_the_per_argument_limit():
    """MAX_ARG_STRLEN bounds each argument, not the command line.

    This is why the budget need not leave room for --model and friends; the
    whole command line is bounded separately by ARG_MAX, which is far larger.
    """
    big = "x" * (MAX_ARG_STRLEN - 1)
    subprocess.run(
        ["/bin/true", big, "--model", "gemini-3.8-flash", "--print-timeout", "600s"],
        check=True,
    )


# --- the budget itself ----------------------------------------------------


def test_budget_leaves_a_margin_below_the_kernel_limit():
    assert AGY_ARGV_BUDGET < MAX_ARG_STRLEN
    margin = MAX_ARG_STRLEN - AGY_ARGV_BUDGET
    assert margin >= 1024, "too little room to fail cleanly rather than with E2BIG"


def test_distiller_and_member_prompts_share_one_budget():
    """Two callers, one kernel limit: raising one alone just moves the wall."""
    assert DEFAULT_MAX_PROMPT_BYTES == AGY_ARGV_BUDGET
    assert AGY_ARGV_LIMIT == AGY_ARGV_BUDGET
