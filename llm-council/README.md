# LLM Council

Karpathy's [llm-council](https://github.com/karpathy/llm-council) rebuilt on subscription CLIs
instead of metered OpenRouter calls, hosted in a [Herdr](https://herdr.dev) multi-pane session.

Three members answer the same question, peer-review each other's answers anonymously, and a
chairman synthesises the result.

| Role | Command | Model |
| --- | --- | --- |
| Member A | `claude -p` | `claude-opus-5-5` |
| Member B | `codex exec` | `gpt-6-astra` |
| Member C | `agy -p` | `gemini-4-argon` |
| Chairman | `claude -p` | `claude-fable-5-1` |
| Distiller | `agy -p` | `gemini-3.8-flash` |

Each council command sets its model explicitly. Member A uses Opus 5.5. The chairman uses
Fable 5.1, the most capable model, because synthesis reads every answer and every review.
Both need a Claude Code build that recognises these model IDs; see the
[model configuration docs](https://code.claude.com/docs/en/model-config).

`agy` rejects a bare model name, so both agy calls also pass `--effort`: member C reasons
at `high`, matching codex; the distiller runs at `low` because it only compresses.

## Status

- [x] 1. Firecrawl CLI installed and authenticated
- [x] 2. Herdr pane API verified
- [x] 3. Shared retrieval + distillation (`council/research.py`)
- [x] 4. Pane scaffolding (`council/panes.py`)
- [x] 5. Stage 1 fan-out with sentinel files (`council/fanout.py`)
- [x] 6. Stages 2–3, shuffle, chairman (`council/deliberate.py`)
- [x] 7. End-to-end run (`council/run.py`)

## Usage

```bash
council run "your question"                      # no retrieval
council run "your question" --web --limit 5      # shared evidence brief first
council run "recent news" --web --tbs qdr:w      # narrow the window to a week
council run "$(cat prompt.xml)" --web --query "short searchable terms"
council run "history of X" --web --tbs none      # search all time, not just recent
council run "your question" --debate             # advocate / critic / neutral stances
council run "your question" --seed 42            # reproducible anonymisation shuffle
council run "your question" --close-panes        # tear the layout down when finished
```

A run without retrieval takes roughly 90 seconds end to end.

## Stage 0: shared retrieval

All members reason over one identical evidence brief, so disagreement reflects interpretation
rather than which documents each happened to retrieve. Full scrapes are kept beside the brief
so a member can consult a complete source when the brief looks thin.

```bash
council research "your question" --limit 5 --run-dir /tmp/council-runs/demo
council research "recent news" --limit 5 --tbs qdr:w --run-dir /tmp/council-runs/demo
```

## The vault note

Every run also writes the synthesis as a note in the personal Obsidian vault, at
`/config/Documents/Harvinder remote vault/copilot/LLM Council/` (`VAULT_DIR` in
`council/vault.py`). Only the note goes there. Prompts, answers, reviews and JSON stay in the run
directory, which the note links back to.

- **Name:** `YYYY-MM-DD HHMM <question>.md`, built from the question's first line with characters
  Obsidian rejects removed and capped at 80 characters. An existing note is never overwritten; a
  clash gets a ` 2`, ` 3` suffix.
- **Content:** YAML frontmatter, then `final.md` unchanged, then a `Run directory` line.
- **Frontmatter:** `created`, `tags: [llm-council]`, `type: LLM Council`, `question`, `panel`,
  `chair`, `ranking` (member names, best first; `[]` when no ballot was usable, since the
  alphabetical fallback is not a verdict), `mean_similarity` and `run_dir`.

The `LLM Council` folder is created if absent, but the vault is not: if the vault is not mounted
the note fails rather than building a stray directory tree. Like a failed export, a failed note is
reported and does not fail the run.

## Exporting the synthesis

Every run already writes two markdown files: `<run-dir>/stage3/chairman.md` is the bare
synthesis, `<run-dir>/final.md` adds the provenance footer. `--export` copies the latter
somewhere useful:

```bash
council run "question" --export ~/Documents/synthesis.md      # copied
council run "question" --export ~/Documents/synthesis.docx    # converted
```

Only `.md`, `.markdown` and `.docx` are supported, validated **before** the run starts so a bad
path fails immediately rather than after several minutes of work. Parent directories are created.

Conversion runs `uv run --with pypandoc_binary`, which fetches a wheel carrying the pandoc binary
into a temporary environment and caches it. **The project itself declares no dependency on
pandoc** — the cost is a ~33 MB download the first time a docx is produced. LibreOffice is not
used: it has no markdown import filter and refuses to load a `.md` file at all.

A failed conversion is reported but does not fail the run — `final.md` is intact either way, and
throwing away a completed council over a format conversion would be absurd.

## Attachments

```bash
council run "does this design hold up?" --attach spec.docx --attach notes.md
council run "summarise the risks" --web --attach report.pdf
```

Supported: `.md`, `.markdown`, `.txt`, `.pdf`, `.docx`, `.doc`, `.odt`, `.rtf`, `.xlsx`, `.xls`,
`.html`, `.htm`. Repeatable.

Markdown and text are read straight off disk. Everything else goes through `firecrawl parse`,
which returns clean markdown and costs credits. Converted text is kept in
`<run-dir>/attachments/A1-<slug>.md`.

Attachments are tagged `[A1]`, `[A2]` and inlined into every member's prompt, so access stays
symmetric — no member needs file tools and none has an advantage over another. They are repeated
in the stage 2 review prompt, because a reviewer that cannot see an attached document cannot
check whether a claim about it is accurate.

Paths are validated before retrieval starts, so a typo fails immediately rather than after
credits and quota have been spent. A document that yields no text — a scanned PDF, say — is an
error rather than a silently empty attachment.

Each document is capped at `DEFAULT_MAX_CHARS_PER_ATTACHMENT` (20 000). Gemini's whole prompt is
still bounded by `MAX_ARG_STRLEN`, so several large attachments plus a brief will raise
`PromptTooLargeError` at dispatch rather than silently truncating.

Images are deliberately unsupported: only Codex can attach them natively, so an image would be
visible to one member and invisible to the other two.

### Recency window

Retrieval is filtered to the **last 180 days by default**. The members are themselves trained on
historical content, so the marginal value of a live search is what postdates their cutoffs;
unfiltered search mostly returns older, better-linked pages the models already know.

The window is `DEFAULT_RECENCY_DAYS` in `research.py`, expressed as an explicit `tbs` date range
(`cdr:1,cd_min:…,cd_max:…`) rather than `qdr:m6`, so it is exactly 180 days rather than
"about six months". Override with any `tbs` value, or pass `none`/`off`/`any` to search all time.

Note that **Firecrawl silently ignores an unrecognised `tbs`** rather than erroring — a typo
degrades to an unfiltered search with no warning. Verified directly: `qdr:zzz9` returns results
as though no filter were set.

**Firecrawl also writes no output file when a search matches nothing**, printing `No results
found.` and exiting 0. A zero exit therefore does not imply the file exists, and a missing file is
read as an empty result set rather than as a failure to run.

### Separating the search query from the question

`--web` searches on the question itself, which fails once the question is a long structured
prompt: search backends cap query length well below a full prompt, so a multi-paragraph or
XML-style question reliably matches nothing. `--query` gives retrieval its own short query while
the panel still receives the full prompt:

```bash
council run "$(cat long-prompt.xml)" --web --query "python-pptx slide master branding"
```

The distiller still sees the real question, so it knows what to keep from the sources; only the
search string changes. `--query` without `--web` is rejected rather than ignored, since nothing
would retrieve it.

### Escalating retrieval

`--limit` is a **target of sources that actually scraped**, not a number of results to request.
Hosts that block scraping — Reddit and LinkedIn among them — return an empty body and are
discarded, so a single search routinely yields fewer usable sources than asked for.

Retrieval therefore retries at increasing depth (5, 10, 20 … capped at Firecrawl's maximum of
100), keeping only URLs not already seen, until the target is met or `--max-attempts` (default
15) is spent. Retrying at the *same* depth would be pointless: Firecrawl returns a
deterministic ranked list, so only a deeper request surfaces new URLs.

It stops early when at maximum depth with nothing new, since further attempts would only spend
credits. If the target is not reached, the run continues on what it found and warns.

**This costs materially more than a single search.** A live query needing 2 attempts (depths 5
then 10) used 15 credits to reach 5 sources, against 3 credits for the single shallow search
that had previously yielded 1 source.

Output:

```
<run-dir>/research/
  search.json                # raw Firecrawl payload
  raw/S1-<slug>.md ...       # full untruncated scrapes with front matter
  brief.md                   # distilled, [S1]-cited, canonical source list
```

Source URLs in `brief.md` are rendered from the Firecrawl payload, never from distiller output,
so a member cannot be pointed at a hallucinated link.

## Pane scaffolding

```bash
council panes create --cwd <run-dir> --state <run-dir>/panes.json --replace
council panes close
```

Layout, one pane per role:

```
+-----------------------------------------+
|              orchestrator               |
+-------------+-------------+-------------+
|  A · opus   |  B · codex  |  C · gemini |
+-------------+-------------+-------------+
|          chair · claude fable           |
+-----------------------------------------+
```

Herdr's `--ratio` is the fraction retained by the *source* pane; the new pane takes the
remainder. `--replace` closes an existing tab of the same label first, so re-running never
leaves duplicates. `--state` writes the role-to-pane map for later stages to re-read.

Panes are created with `--no-focus`: a council run must not yank you out of the tab you are
working in.

## Stage 1: fan-out

All three members receive the identical prompt and run in parallel, each in its own pane.
A live run of three members completed in about 20 seconds.

Each member's pane runs a generated shell script that:

1. prints a banner naming the member,
2. runs the agent headless, piping output through `tee` so the pane shows the work *and* a
   clean copy lands on disk,
3. writes the exit status to `<role>.status`,
4. touches `<role>.done`.

The orchestrator polls for the sentinels. It never reads rendered terminal output, so ANSI
codes, spinners and line wrapping cannot corrupt an answer.

```
<run-dir>/stage1/
  <role>.prompt.md   <role>.sh      <role>.md
  <role>.status      <role>.done    <role>.log
```

Two properties the scripts must preserve, both covered by tests:

- **No `set -e`.** An aborting script would never reach the sentinel and the orchestrator
  would block until the stage timed out.
- **`set -o pipefail`.** The agent's exit status has to survive the pipe into `tee`.

Stale sentinels are deleted at dispatch, so a re-run can never collect the previous answer.

### Brief versus recollection

Members cannot retrieve anything — Claude runs `--allowedTools ""`, Codex a read-only sandbox — so
the brief is the whole of their *fresh* information. Their training data is unavoidably active
anyway, supplying reasoning, framing and judgement. Two rules govern the interaction, and they
pull in opposite directions deliberately:

1. **On contradiction, the brief wins.** It was retrieved inside the recency window; training data
   is the older source. Members are asked to say so explicitly and name what they previously
   believed, since a contradiction is useful signal rather than something to smooth over.
2. **Where the brief is silent or thin, use prior knowledge.** Five scraped pages are not a survey
   and coverage may be poor or skewed. Such claims are marked `[prior]` so a reader can separate
   retrieved fact from recollection.

**The review prompt is coupled to this and must stay in sync.** It previously said "penalise claims
the brief does not support", which would punish exactly the supplementation rule 2 asks for. It now
tells reviewers not to penalise a *marked* claim for lacking brief support, while still penalising
unmarked unsupported claims and unacknowledged contradictions. Changing one prompt without the
other puts the two stages in conflict.

Measured on an A/B against the previous wording, same brief and same members:

| Case | Old wording | New wording |
| --- | --- | --- |
| Real brief, niche topic | states the version correctly, invents nothing | unchanged, plus an explicit "where this contradicts what I knew" section |
| Deliberately thin brief, well-known topic | **named the gaps and answered nothing** ("I'm not going to fill those gaps from memory here") | answered, with claims marked `[prior]` |

Costs: answers run 2–3× longer, and marker discipline varies a lot between members (0 to 17 marks
on comparable answers). A missing `[prior]` marker is therefore not evidence that an answer is
fully grounded.

### Failure modes that are treated as errors

- non-zero exit, reported with the tail of the member's log
- exit 0 with an empty answer, which is how `agy` reports an auto-denied tool permission
- a member that never signals, surfaced as a timeout naming who stalled

## Stages 2 and 3: review and synthesis

Answers are relabelled `A`/`B`/`C` under a fresh shuffle each run before any member sees
them, so a reviewer cannot favour a vendor it recognises. Reviewers also receive the shared
evidence brief — that is what makes peer review meaningful, since a reviewer can check a claim
against the same sources it was drawn from instead of judging on plausibility alone.

Each reviewer ranks all three answers, including its own, and must end with a parseable line:

```
RANKING: B > A > C
```

Ballots are aggregated by Borda count. A reviewer that produces no usable ranking is recorded
in `rankings.json` under `failures` and simply excluded — one bad ballot does not sink the
stage. The chair receives the answers, the aggregate order and the individual ballots, and is
explicitly told the ranking is evidence rather than an instruction to copy the winner.

The label-to-vendor mapping is disclosed in `final.md`. Anonymity exists to stop members
favouring each other, not to keep you in the dark.

```
<run-dir>/
  query.txt  labels.json  rankings.json  divergence.json  panes.json
  research/  stage1/  stage2/  stage3/  final.md
```

### Divergence logging

Every run records the mean pairwise similarity of the stage 1 answers (Jaccard over word
bigrams) to `divergence.json`, so convergence is observable over time rather than guessed at.

**The absolute numbers run low and the 0.75 `converged` threshold is not yet calibrated.**
Observed means on real runs were 0.019–0.095, because two answers must phrase things alike,
not merely discuss the same topic, to score highly. Treat the metric as a relative trend
across runs; the threshold needs retuning against a corpus of real runs before it will fire
usefully.

### Debate mode

`--debate` gives each member a different stance — advocate, critic, neutral — fixed per pane
so the layout stays legible. This deliberately trades the "which model is best" comparison for
coverage of the argument space, so it is an exploratory mode rather than a measurement one.

## Escape hatch (`--escape-hatch`, opt-in)

Lets a member read the full text of a source when the brief is too compressed to settle
something. Requires `--web`, since without retrieval there are no sources to read.

Antigravity's permission matcher accepts **only literal paths** — every glob form was tested
and rejected, including `**`, `*` and bare tool names. A malformed rule silently discards the
entire allow list (the log shows `permissions=<nil>`), so rules must be exact.

Each run therefore mirrors its raw scrapes to fixed slot filenames that a static rule set can
cover:

```
/config/.llm-council/current/raw/S1.md … S10.md
```

Slot `Sn` corresponds to source tag `[Sn]` in the brief. Unused slots are overwritten with a
placeholder each run, so a shorter run cannot leave the previous run's sources sitting at a
readable path where a member might cite them as current.

Enable by adding literal rules to `/config/.gemini/antigravity-cli/settings.json`:

```json
"permissions": {
  "allow": [
    "read_file(/config/.llm-council/current/raw/S1.md)",
    "read_file(/config/.llm-council/current/raw/S2.md)"
  ]
}
```

Grants per member:

| Member | Mechanism | Can use hatch |
| --- | --- | --- |
| A, claude | `--allowedTools "Read"` — never Edit, Write or Bash | yes |
| C, agy | literal `read_file(...)` rules in antigravity settings | yes |
| B, codex | none | **no** |

**Codex is excluded deliberately.** Its filesystem sandbox fails closed in this container — it
cannot read a file even in its own working directory, and neither `-C` nor
`sandbox_permissions=["disk-full-read-access"]` changes that. The only mechanism that would
work is `--sandbox danger-full-access`, which removes sandboxing from every council run in
exchange for an occasional source read. Not worth it, so Codex answers from the brief alone.

That asymmetry is real and is accepted knowingly: two members can consult fuller evidence than
the third, so **peer ranking partly reflects tooling rather than reasoning** when the hatch is
on. To keep it visible, members cite hatch-derived claims as `[Sn:full]` and usage is recorded
to `hatch.json` per run. That marker is self-reported and undercounts a member that forgets it;
it exists to expose the asymmetry, not to audit it. Use plain `council run` without
`--escape-hatch` when you want a clean model comparison.

Known limitations:

- **One run at a time.** The mirror is a single shared location, so concurrent runs clobber it.
- **Sources beyond slot 10** are archived in the run directory but not reachable via the hatch.
- **Members must be told exact paths.** A member that greps for these files is denied; the
  prompt says so explicitly.

## Not yet implemented

**`discoveries/` pooling.** The design called for anything a member retrieves independently to
be pooled so all reviewers see it before stage 2. Members currently have no web search — Claude
runs `--allowedTools ""` and Codex `--sandbox read-only` — so there is nothing to pool. This
becomes relevant only if per-member search is enabled later.

## Tests

```bash
uv run --with pytest pytest tests/ -q
```

Subprocess boundaries take an injectable `runner`, so the suite exercises command construction
and failure handling without spending credits or quota.

`run_council` takes the same seam, and it is the only one it needs: every child process — Herdr,
Firecrawl, `agy`, `uv` — is started through that one callable. A fake that routes on the
executable name and, when Herdr is told to run a member's script, writes the answer and sentinel
that script would have written, drives a whole run in milliseconds. Stage sequencing, the
no-usable-ballot fallback, member failure, stage timeout and the export-failure path are all
covered that way rather than by real runs.

`mirror_dir` is injectable for the same reason: the escape-hatch mirror is a fixed absolute path,
so without it a test of `research()` would overwrite the live mirror.

The vault is protected the same way, but automatically: an autouse fixture in `tests/conftest.py`
points `council.vault.VAULT_DIR` at a temporary folder, so no test can write into the real vault
even if it forgets to pass `vault_dir`.

## Environment constraints

These are verified facts about this machine, not assumptions:

- **`agy` ignores stdin.** Prompts must be passed via argv, so they are bounded by
  `MAX_ARG_STRLEN` — `32 × PAGE_SIZE`, hence 131072 bytes on x86-64. Measured here: 131071 bytes
  execs, 131072 raises `E2BIG`. Distiller and member prompts share one `AGY_ARGV_BUDGET` of
  128000, raising `CorpusTooLargeError` and `PromptTooLargeError` respectively past that.

  The limit is **per argument, not per command line**: `--model` and the rest cost nothing against
  it, and the command line as a whole is bounded separately by `ARG_MAX` (2 MB here), which
  nothing in this project approaches. The 3072-byte margin therefore is not room for other
  arguments — it exists so an oversized prompt fails as a clear error rather than as an opaque
  exec failure inside a pane script. `tests/test_argv_limits.py` pins all of this by actually
  exec'ing at the boundary.
- **`agy` headless denies tools that need approval.** This includes shell commands and file reads.
  A peer review failed when Gemini requested `cal 9 2026` to check a date.
  All council prompts now prohibit shell commands and other tools, except the explicit source reads for `--escape-hatch`.
  Review prompts also tell members to rank the responses without executing the original request.
  Members must state uncertainty when they cannot verify a detail from the supplied text or their own knowledge.

  These prompt instructions guide the model. They do not disable tools at the CLI level.
  Council still rejects empty answers if a model requests a tool that needs approval.
  See the [Antigravity headless permissions documentation](https://antigravity.google/docs/cli/headless/#permissions-in-headless-mode).

  The `raw/` escape hatch still requires literal `read_file(...)` rules under `permissions.allow`
  in `/config/.gemini/antigravity-cli/settings.json`.
- **`agy` exits 0 with empty stdout** when a permission is denied, so empty output is treated as
  an error rather than an empty answer.
- **Firecrawl concurrency is 2**, not the 100 shown in the skill docs.
- **`codex exec` needs `--skip-git-repo-check`** outside a git repository.
- npm's global prefix was repointed to `/config/.local` so `firecrawl` installs without root.
  Reversible with `npm config delete prefix`.

## Cost

Retrieval runs about 5 credits for 3 results and 7–8 for 5. A full council query is 8 agent
invocations: 1 distil + 3 members + 3 reviews + 1 chairman.
