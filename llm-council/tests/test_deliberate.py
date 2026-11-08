# ABOUTME: Tests for stage 2 peer review and stage 3 synthesis.
# ABOUTME: Covers anonymisation, ranking parsing, Borda aggregation and divergence scoring.

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from council.deliberate import (  # noqa: E402
    MalformedRankingError,
    aggregate_rankings,
    assign_labels,
    build_chair_prompt,
    build_review_prompt,
    divergence_report,
    parse_ranking,
    similarity,
)
from council.panes import MEMBER_A, MEMBER_B, MEMBER_C  # noqa: E402

ROLES = [MEMBER_A, MEMBER_B, MEMBER_C]
ANSWERS = {
    MEMBER_A: "The sky is blue because of Rayleigh scattering [S1].",
    MEMBER_B: "Short wavelengths scatter more in the atmosphere [S2].",
    MEMBER_C: "Blue light scatters most, so the sky looks blue [S1].",
}


# --- anonymisation --------------------------------------------------------


def test_assign_labels_is_a_bijection_onto_letters():
    labels = assign_labels(ROLES, rng=random.Random(1))
    assert set(labels) == set(ROLES)
    assert sorted(labels.values()) == ["A", "B", "C"]


def test_assign_labels_is_deterministic_for_a_given_seed():
    assert assign_labels(ROLES, rng=random.Random(7)) == assign_labels(ROLES, rng=random.Random(7))


def test_assign_labels_actually_shuffles_across_runs():
    """A fixed mapping would let a model learn which letter is which vendor."""
    seen = {tuple(sorted(assign_labels(ROLES, rng=random.Random(s)).items())) for s in range(40)}
    assert len(seen) > 1


# --- review prompt --------------------------------------------------------


def test_review_prompt_presents_every_answer_under_its_label():
    labels = {MEMBER_A: "B", MEMBER_B: "C", MEMBER_C: "A"}
    prompt = build_review_prompt("why?", ANSWERS, labels, brief=None)
    for role, label in labels.items():
        assert f"Response {label}" in prompt
        assert ANSWERS[role] in prompt


def test_review_prompt_never_leaks_role_or_vendor_identity():
    """Anonymity is the point: a reviewer must not know who wrote what."""
    labels = {MEMBER_A: "B", MEMBER_B: "C", MEMBER_C: "A"}
    prompt = build_review_prompt("why?", ANSWERS, labels, brief=None).lower()
    for leak in ("member_a", "member_b", "member_c", "claude", "codex", "gemini", "opus", "gpt"):
        assert leak not in prompt


def test_review_prompt_orders_responses_by_label_not_by_role():
    labels = {MEMBER_A: "C", MEMBER_B: "A", MEMBER_C: "B"}
    prompt = build_review_prompt("why?", ANSWERS, labels, brief=None)
    assert prompt.index("Response A") < prompt.index("Response B") < prompt.index("Response C")


def test_review_prompt_includes_the_shared_brief_so_claims_can_be_checked():
    prompt = build_review_prompt("why?", ANSWERS, {r: l for r, l in zip(ROLES, "ABC")}, brief="[S1] fact")
    assert "[S1] fact" in prompt


def test_review_prompt_does_not_penalise_declared_prior_knowledge():
    """Members are told to supplement thin coverage; reviewers must not punish them for it."""
    prompt = build_review_prompt(
        "why?", ANSWERS, {r: l for r, l in zip(ROLES, "ABC")}, brief="[S1] fact"
    )
    assert "[prior]" in prompt


def test_review_prompt_still_penalises_undeclared_unsupported_claims():
    prompt = build_review_prompt(
        "why?", ANSWERS, {r: l for r, l in zip(ROLES, "ABC")}, brief="[S1] fact"
    ).lower()
    assert "unmarked" in prompt or "undeclared" in prompt or "without marking" in prompt


def test_review_prompt_demands_a_parseable_ranking_line():
    prompt = build_review_prompt("why?", ANSWERS, {r: l for r, l in zip(ROLES, "ABC")}, brief=None)
    assert "RANKING:" in prompt


# --- ranking parsing ------------------------------------------------------


@pytest.mark.parametrize(
    "text",
    [
        "RANKING: B > A > C",
        "blah blah\n\nRANKING: B > A > C",
        "ranking: b > a > c",
        "**RANKING: B > A > C**",
        "RANKING: B>A>C",
        "RANKING:  B  >  A  >  C  .",
    ],
)
def test_parse_ranking_accepts_realistic_formatting(text):
    assert parse_ranking(text, expected={"A", "B", "C"}) == ["B", "A", "C"]


def test_parse_ranking_uses_the_last_ranking_line():
    """Models sometimes restate a draft ranking before committing to a final one."""
    text = "RANKING: A > B > C\nOn reflection...\nRANKING: C > B > A"
    assert parse_ranking(text, expected={"A", "B", "C"}) == ["C", "B", "A"]


def test_parse_ranking_rejects_a_missing_line():
    with pytest.raises(MalformedRankingError):
        parse_ranking("I prefer B, then A, then C.", expected={"A", "B", "C"})


def test_parse_ranking_rejects_an_incomplete_ranking():
    with pytest.raises(MalformedRankingError):
        parse_ranking("RANKING: B > A", expected={"A", "B", "C"})


def test_parse_ranking_rejects_duplicates():
    with pytest.raises(MalformedRankingError):
        parse_ranking("RANKING: B > B > A", expected={"A", "B", "C"})


def test_parse_ranking_rejects_unknown_labels():
    with pytest.raises(MalformedRankingError):
        parse_ranking("RANKING: B > A > D", expected={"A", "B", "C"})


# --- aggregation ----------------------------------------------------------


def test_aggregate_rankings_uses_borda_counts():
    rankings = {
        "A": ["B", "A", "C"],
        "B": ["B", "C", "A"],
        "C": ["B", "A", "C"],
    }
    order, scores = aggregate_rankings(rankings)
    assert order[0] == "B"
    assert scores["B"] == 6  # first place on all three ballots, 2 points each


def test_aggregate_rankings_breaks_ties_deterministically():
    rankings = {"A": ["A", "B", "C"], "B": ["B", "A", "C"]}
    first, _ = aggregate_rankings(rankings)
    second, _ = aggregate_rankings(rankings)
    assert first == second


def test_aggregate_rankings_tolerates_a_reviewer_that_failed():
    """One unusable ballot must not sink the whole stage."""
    order, scores = aggregate_rankings({"A": ["C", "B", "A"]})
    assert order[0] == "C"


# --- chair prompt ---------------------------------------------------------


def test_chair_prompt_carries_answers_and_the_peer_verdict():
    labels = {r: l for r, l in zip(ROLES, "ABC")}
    prompt = build_chair_prompt(
        "why?", ANSWERS, labels, order=["B", "A", "C"], brief=None, rankings={"A": ["B", "A", "C"]}
    )
    assert "why?" in prompt
    for answer in ANSWERS.values():
        assert answer in prompt
    assert "B > A > C" in prompt


def test_chair_prompt_stays_anonymous():
    labels = {r: l for r, l in zip(ROLES, "ABC")}
    prompt = build_chair_prompt(
        "why?", ANSWERS, labels, order=["A", "B", "C"], brief=None, rankings={}
    ).lower()
    for leak in ("claude", "codex", "gemini", "member_a"):
        assert leak not in prompt


def test_chair_prompt_tells_the_chair_to_judge_not_to_defer():
    labels = {r: l for r, l in zip(ROLES, "ABC")}
    prompt = build_chair_prompt(
        "why?", ANSWERS, labels, order=["A", "B", "C"], brief=None, rankings={}
    )
    assert "ranking" in prompt.lower()
    # The peer ranking is evidence, not an instruction to copy the winner.
    assert "not" in prompt.lower()


# --- divergence -----------------------------------------------------------


def test_similarity_of_identical_text_is_one():
    assert similarity("the cat sat on the mat", "the cat sat on the mat") == 1.0


def test_similarity_of_disjoint_text_is_zero():
    assert similarity("alpha beta gamma", "delta epsilon zeta") == 0.0


def test_similarity_is_symmetric():
    a, b = "the cat sat down", "the cat stood up"
    assert similarity(a, b) == similarity(b, a)


def test_similarity_ignores_case_and_punctuation():
    assert similarity("The Cat, sat!", "the cat sat") == 1.0


def test_divergence_report_scores_every_pair():
    report = divergence_report(ANSWERS, {r: l for r, l in zip(ROLES, "ABC")})
    assert len(report["pairs"]) == 3
    assert set(report["pairs"][0]) >= {"a", "b", "similarity"}


def test_divergence_report_flags_high_convergence():
    """If members agree almost verbatim, the council is adding little."""
    same = {role: "identical text here" for role in ROLES}
    report = divergence_report(same, {r: l for r, l in zip(ROLES, "ABC")})
    assert report["mean_similarity"] == 1.0
    assert report["converged"] is True
