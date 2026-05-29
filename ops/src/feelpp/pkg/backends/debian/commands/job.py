from __future__ import annotations

import argparse
import json

from ....workspace import ensure_workspace, write_job_manifest
from .common import add_context_arguments, add_runtime_arguments, context_from_args


def command_job_init(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    ensure_workspace(context)
    write_job_manifest(context, state="initialized")
    print(
        json.dumps(
            {"context": context.as_dict(), "job_manifest": str(context.job_manifest_path)},
            indent=2,
        )
    )
    return 0


def register_job_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    job_parser = subparsers.add_parser("job", help="Job workspace commands")
    job_subparsers = job_parser.add_subparsers(dest="job_command", required=True)
    job_init_parser = job_subparsers.add_parser("init", help="Initialize the packaging workspace")
    add_context_arguments(job_init_parser)
    add_runtime_arguments(job_init_parser, default_engine="host")
    job_init_parser.set_defaults(func=command_job_init)
