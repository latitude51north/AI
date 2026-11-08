# ABOUTME: Stage 2 peer review and stage 3 synthesis for the LLM Council.
# ABOUTME: Handles anonymisation, ranking parsing, Borda aggregation and divergence scoring.

"""Deliberation.

Stage 1 answers are relabelled A/B/C under a per-run shuffle before any member
sees them, so a reviewer cannot favour a vendor it recognises. Reviewers also
receive the shared evidence brief, which is what makes peer review meaningful:
a reviewer can check a claim against the same sources the claim was drawn from,
rather than judging on plausibility alone.
"""

from __future__ import annotations

import random
import re
from typing import Iterable, Mapping, Sequence

from council.prompts import headless_instructions

# A reviewer must end with a line in this form so the ballot can be parsed.
RANKING_PATTERN = re.compile(r"RANKING\s*:\s*([A-Za-z](?:\s*>\s*[A-Za-z])+)", re.IGNORECASE)

# Above this mean pairwise similarity the members are saying the same thing and
# the council is adding little over a single call.
CONVERGENCE_THRESHOLD = 0.75

_WORD = re.compile(r"[a-z0-9]+")


class MalformedRankingError(RuntimeError):
    """A reviewer did not produce a usable ranking line."""


# --- anonymisation --------------------------------------------------------


def assign_labels(roles: Sequence[str], rng: random.Random | None = None) -> dict[str, str]:
    """Map each role to a letter under a fresh per-run shuffle."""
    rng = rng or random.Random()
    letters = [chr(ord("A") + i) for i in range(len(roles))]
    rng.shuffle(letters)
    return dict(zip(roles, letters))


def _by_label(answers: Mapping[str, str], labels: Mapping[str, str]) -> list[tuple[str, str]]:
    """Answers ordered by label, so presentation order carries no role information."""
    return sorted(((labels[role], text) for role, text in answers.items()), key=lambda p: p[0])


# --- stage 2 --------------------------------------------------------------


def build_review_prompt(
    question: str,
    answers: Mapping[str, str],
    labels: Mapping[str, str],
    brief: str | None,
    attachments: str | None = None,
) -> str:
    blocks = [
        f"### Response {label}\n\n{text}" for label, text in _by_label(answers, labels)
    ]
    letters = [label for label, _ in _by_label(answers, labels)]
    example = " > ".join(letters)

    sections = [
        "Three responses to the same question are shown below. They were written "
        "independently and are presented anonymously.",
        "",
        headless_instructions(),
        "Review and rank the responses. Do not carry out the original request yourself.",
        "",
        f"QUESTION:\n{question}",
    ]
    if attachments:
        # Reviewers need the attachments too, or they cannot check whether a
        # claim about an attached document is accurate.
        sections += [
            "",
            "All three authors were given the attached documents below. Use them to "
            "check whether a response's claims about those documents are accurate.",
            "",
            f"ATTACHED DOCUMENTS:\n{attachments}",
        ]
    if brief:
        sections += [
            "",
            "All three authors worked from the shared evidence brief below. Use it to "
            "check whether a response's claims are actually supported.",
            "",
            f"EVIDENCE BRIEF:\n{brief}",
        ]
    sections += [
        "",
        "\n\n".join(blocks),
        "",
        "Rank the responses on accuracy first and insight second. Reward a response "
        "that is correctly cautious over one that is confidently wrong. Ignore length, "
        "formatting and tone.",
        "",
        "On claims that go beyond the brief: the authors were told to use their own "
        "knowledge where the brief is thin, and to mark such claims [prior]. Do not "
        "penalise a marked claim merely for being unsupported by the brief — judge "
        "whether it is correct and whether it earns its place. Do penalise an unmarked "
        "claim that the brief does not support, and penalise any claim that contradicts "
        "the brief without acknowledging the contradiction.",
        "",
        "Give a short justification for each placement, then end your reply with a single "
        f"line in exactly this form:\n\nRANKING: {example}",
    ]
    return "\n".join(sections) + "\n"


def parse_ranking(text: str, expected: Iterable[str]) -> list[str]:
    """Pull the ballot out of a reviewer's reply.

    The last matching line wins: a reviewer that reconsiders mid-reply should be
    taken at its final word.
    """
    expected = {label.upper() for label in expected}
    matches = RANKING_PATTERN.findall(text or "")
    if not matches:
        raise MalformedRankingError("no RANKING line found in the reviewer's reply")

    order = [part.strip().upper() for part in re.split(r">", matches[-1])]

    if len(order) != len(expected):
        raise MalformedRankingError(f"ranking {order} does not cover all of {sorted(expected)}")
    if len(set(order)) != len(order):
        raise MalformedRankingError(f"ranking {order} repeats a label")
    if set(order) != expected:
        raise MalformedRankingError(f"ranking {order} does not match {sorted(expected)}")
    return order


def aggregate_rankings(rankings: Mapping[str, Sequence[str]]) -> tuple[list[str], dict[str, int]]:
    """Borda count across whatever ballots parsed successfully.

    A reviewer whose ballot could not be parsed is simply absent; one bad ballot
    should not invalidate the stage.
    """
    scores: dict[str, int] = {}
    for ballot in rankings.values():
        size = len(ballot)
        for position, label in enumerate(ballot):
            scores[label] = scores.get(label, 0) + (size - 1 - position)

    # Sort by score descending, then label ascending so ties are stable.
    order = sorted(scores, key=lambda label: (-scores[label], label))
    return order, scores


# --- stage 3 --------------------------------------------------------------


def build_chair_prompt(
    question: str,
    answers: Mapping[str, str],
    labels: Mapping[str, str],
    order: Sequence[str],
    brief: str | None,
    rankings: Mapping[str, Sequence[str]],
) -> str:
    blocks = [f"### Response {label}\n\n{text}" for label, text in _by_label(answers, labels)]
    ballots = "\n".join(f"- a reviewer ranked: {' > '.join(b)}" for b in rankings.values())

    sections = [
        "You are the chair of a panel. Three anonymous responses to the same question "
        "are below, along with how the panel ranked them.",
        "",
        headless_instructions(),
        "",
        f"QUESTION:\n{question}",
    ]
    if brief:
        sections += ["", f"EVIDENCE BRIEF:\n{brief}"]
    sections += [
        "",
        "\n\n".join(blocks),
        "",
        "PANEL RANKING (best first): " + " > ".join(order),
    ]
    if ballots:
        sections += ["", "Individual ballots:", ballots]
    sections += [
        "",
        "Write the single best answer to the question. Take the strongest reasoning from "
        "wherever it appears, including from a response the panel ranked last. The panel "
        "ranking is evidence about quality, not an instruction to copy the winner: if the "
        "panel favoured a response you judge to be wrong, say so and explain why.",
        "",
        "Where the responses genuinely disagree, resolve it if the evidence allows and "
        "state the disagreement plainly if it does not. Do not mention the panel, the "
        "responses, or the ranking in your answer. Write it as a direct answer to the "
        "question with no preamble.",
    ]
    return "\n".join(sections) + "\n"


# --- divergence -----------------------------------------------------------


def _shingles(text: str) -> set[tuple[str, str]]:
    words = _WORD.findall((text or "").lower())
    return set(zip(words, words[1:]))


def similarity(a: str, b: str) -> float:
    """Jaccard overlap of word bigrams.

    Bigrams rather than single words so that shared vocabulary on a shared topic
    does not read as agreement; two answers must phrase things alike to score high.
    """
    sa, sb = _shingles(a), _shingles(b)
    if not sa and not sb:
        return 1.0
    if not sa or not sb:
        return 0.0
    return len(sa & sb) / len(sa | sb)


def divergence_report(answers: Mapping[str, str], labels: Mapping[str, str]) -> dict:
    """Score how much the members actually differed.

    Logged every run so convergence is something you can observe over time rather
    than guess at.
    """
    labelled = _by_label(answers, labels)
    pairs = []
    for i, (label_a, text_a) in enumerate(labelled):
        for label_b, text_b in labelled[i + 1 :]:
            pairs.append(
                {"a": label_a, "b": label_b, "similarity": round(similarity(text_a, text_b), 4)}
            )

    mean = sum(p["similarity"] for p in pairs) / len(pairs) if pairs else 0.0
    return {
        "pairs": pairs,
        "mean_similarity": round(mean, 4),
        "converged": mean >= CONVERGENCE_THRESHOLD,
        "threshold": CONVERGENCE_THRESHOLD,
    }
