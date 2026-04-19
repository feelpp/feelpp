from __future__ import annotations

import argparse
import json
from pathlib import Path

from ..localrepo import stage_outputs
from .common import add_context_arguments, add_runtime_arguments, context_from_args


def command_repo_stage(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    result = stage_outputs(context, Path(args.result_dir).expanduser().resolve())
    print(json.dumps(result, indent=2))
    return 0


def register_repo_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    repo_parser = subparsers.add_parser("repo", help="Local repository commands")
    repo_subparsers = repo_parser.add_subparsers(dest="repo_command", required=True)
    repo_stage_parser = repo_subparsers.add_parser(
        "stage",
        help="Stage build outputs into the local repo",
    )
    add_context_arguments(repo_stage_parser)
    add_runtime_arguments(repo_stage_parser, default_engine="host")
    repo_stage_parser.add_argument("result_dir", help="Result directory to stage")
    repo_stage_parser.set_defaults(func=command_repo_stage)
