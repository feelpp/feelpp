from __future__ import annotations

import argparse
import json

from .common import add_context_arguments, add_runtime_arguments, context_from_args, load_plan


def command_inspect_plan(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    plan = load_plan(args, context)
    payload = {
        "context": context.as_dict(),
        "plan": plan.as_dict(),
        "manifest_path": str(context.manifest_path),
    }
    print(json.dumps(payload, indent=2))
    return 0


def register_inspect_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    inspect_parser = subparsers.add_parser("inspect", help="Inspection commands")
    inspect_subparsers = inspect_parser.add_subparsers(dest="inspect_command", required=True)
    inspect_plan_parser = inspect_subparsers.add_parser("plan", help="Show the resolved component plan")
    add_context_arguments(inspect_plan_parser)
    add_runtime_arguments(inspect_plan_parser, default_engine="host")
    inspect_plan_parser.add_argument("--components", help="Comma-separated component list override")
    inspect_plan_parser.add_argument(
        "--skip-component",
        action="append",
        default=[],
        help="Component to skip",
    )
    inspect_plan_parser.add_argument(
        "--skip-text",
        default="",
        help="Raw skip text such as a commit message",
    )
    inspect_plan_parser.add_argument(
        "--publish",
        action="store_true",
        help="Resolve plan with publish enabled",
    )
    inspect_plan_parser.set_defaults(func=command_inspect_plan)
