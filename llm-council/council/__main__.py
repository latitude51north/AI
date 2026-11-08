# ABOUTME: Command-line entrypoint for the LLM Council.
# ABOUTME: Dispatches subcommands to the stage modules.

from __future__ import annotations

import sys
from typing import Sequence

USAGE = """council <command> [options]

Commands:
  run        Ask the full council a question, end to end
  research   Run the shared retrieval stage and write an evidence brief
  panes      Create or tear down the Herdr pane layout

Run `council <command> --help` for command options.
"""


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(USAGE)
        return 0

    command, rest = argv[0], argv[1:]

    if command == "run":
        from council.run import main as run_main

        return run_main(rest)

    if command == "research":
        from council.research import main as research_main

        return research_main(rest)

    if command == "panes":
        from council.panes import main as panes_main

        return panes_main(rest)

    print(f"unknown command: {command}\n\n{USAGE}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
