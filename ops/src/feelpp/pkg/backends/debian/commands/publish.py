from __future__ import annotations

import argparse
from pathlib import Path

from ....publish import publish_cleanup, publish_snapshot
from .common import add_context_arguments, add_runtime_arguments, context_from_args


def command_publish_snapshot(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    input_dir = Path(args.input_dir).expanduser().resolve() if args.input_dir else None
    publish_snapshot(context, input_dir=input_dir, dry_run=args.dry_run)
    return 0


def command_publish_cleanup(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    publish_cleanup(context, dry_run=args.dry_run)
    return 0


def register_publish_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    publish_parser = subparsers.add_parser("publish", help="Publish commands")
    publish_subparsers = publish_parser.add_subparsers(dest="publish_command", required=True)

    publish_snapshot_parser = publish_subparsers.add_parser(
        "snapshot",
        help="Create or switch to a published snapshot",
    )
    add_context_arguments(publish_snapshot_parser)
    add_runtime_arguments(publish_snapshot_parser)
    publish_snapshot_parser.add_argument(
        "--input-dir",
        help="Directory containing package artifacts to publish",
    )
    publish_snapshot_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    publish_snapshot_parser.set_defaults(func=command_publish_snapshot)

    publish_switch_parser = publish_subparsers.add_parser(
        "switch",
        help="Switch the published endpoint using the current snapshot logic",
    )
    add_context_arguments(publish_switch_parser)
    add_runtime_arguments(publish_switch_parser)
    publish_switch_parser.add_argument(
        "--input-dir",
        help="Directory containing package artifacts to publish",
    )
    publish_switch_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    publish_switch_parser.set_defaults(func=command_publish_snapshot)

    publish_cleanup_parser = publish_subparsers.add_parser(
        "cleanup",
        help="Remove prerelease packages whose final release already exists, then republish",
    )
    add_context_arguments(publish_cleanup_parser)
    add_runtime_arguments(publish_cleanup_parser)
    publish_cleanup_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    publish_cleanup_parser.set_defaults(func=command_publish_cleanup)
