from __future__ import annotations

import argparse

from ....build import build_chain, build_component
from ....publish import publish_snapshot
from .common import add_context_arguments, add_runtime_arguments, context_from_args, load_plan


def command_build_component(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    build_component(
        context,
        args.name,
        dry_run=args.dry_run,
        skip_pbuilder_prepare=args.skip_pbuilder_prepare,
    )
    return 0


def command_build_chain(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    plan = load_plan(args, context)
    build_chain(
        context,
        plan,
        dry_run=args.dry_run,
        prepare_pbuilder=not args.skip_pbuilder_prepare,
    )
    if args.publish and plan.publish_enabled:
        publish_snapshot(context, dry_run=args.dry_run)
    return 0


def register_build_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    build_parser = subparsers.add_parser("build", help="Build commands")
    build_subparsers = build_parser.add_subparsers(dest="build_command", required=True)

    build_component_parser = build_subparsers.add_parser("component", help="Build one source component")
    add_context_arguments(build_component_parser)
    add_runtime_arguments(build_component_parser)
    build_component_parser.add_argument("name", help="Component name")
    build_component_parser.add_argument(
        "--skip-pbuilder-prepare",
        action="store_true",
        help="Assume the pbuilder base is already prepared",
    )
    build_component_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    build_component_parser.set_defaults(func=command_build_component)

    build_chain_parser = build_subparsers.add_parser("chain", help="Build the ordered component chain")
    add_context_arguments(build_chain_parser)
    add_runtime_arguments(build_chain_parser)
    build_chain_parser.add_argument("--components", help="Comma-separated component list override")
    build_chain_parser.add_argument(
        "--skip-component",
        action="append",
        default=[],
        help="Component to skip",
    )
    build_chain_parser.add_argument(
        "--skip-text",
        default="",
        help="Raw skip text such as a commit message",
    )
    build_chain_parser.add_argument(
        "--skip-pbuilder-prepare",
        action="store_true",
        help="Assume the pbuilder base is already prepared",
    )
    build_chain_parser.add_argument(
        "--publish",
        action="store_true",
        help="Publish after the build chain completes",
    )
    build_chain_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    build_chain_parser.set_defaults(func=command_build_chain)
