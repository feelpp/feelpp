from __future__ import annotations

import argparse

from ....pbuilder import prepare_base
from .common import add_context_arguments, add_runtime_arguments, context_from_args


def command_pbuilder_prepare(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    prepare_base(context, dry_run=args.dry_run)
    return 0


def register_pbuilder_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    pbuilder_parser = subparsers.add_parser("pbuilder", help="Pbuilder commands")
    pbuilder_subparsers = pbuilder_parser.add_subparsers(dest="pbuilder_command", required=True)
    pbuilder_prepare_parser = pbuilder_subparsers.add_parser(
        "prepare",
        help="Prepare the pbuilder base",
    )
    add_context_arguments(pbuilder_prepare_parser)
    add_runtime_arguments(pbuilder_prepare_parser)
    pbuilder_prepare_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    pbuilder_prepare_parser.set_defaults(func=command_pbuilder_prepare)
