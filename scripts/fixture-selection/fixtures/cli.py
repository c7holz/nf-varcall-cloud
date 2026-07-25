"""Command-line entry point.

Run from ``scripts/fixture-selection/``::

    python3 -m fixtures.cli reproduce

The audit is a single command, ``reproduce``: it re-derives the window and the two
variants from the primary sources and diffs the result against the committed
``results/`` TSVs. There is no separate source-verification step — the sources are
version-pinned and immutable, so re-deriving the numbers already exercises the fetch,
and a source that rots or changes shape makes it *raise* rather than lie.
"""

from __future__ import annotations

import argparse
import sys


def reproduce() -> int:
    """Re-derive the decision and diff it against the committed ``results/``."""
    raise SystemExit(
        "reproduce is not implemented yet (phase 4).\n"
        "The primary-source layer (fixtures.sources) is in place; the site and\n"
        "window filters it feeds, and the results/ TSVs it diffs against, are still\n"
        "to come. See scripts/fixture-selection/README.md for the phase table."
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="fixtures", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser(
        "reproduce",
        help="re-derive the decision from primary sources and diff committed results/",
    )
    args = parser.parse_args(argv)

    if args.command == "reproduce":
        return reproduce()
    parser.error(f"unknown command {args.command}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
