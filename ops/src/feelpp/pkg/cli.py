from __future__ import annotations

import argparse
import os
import subprocess
import sys

from feelpp.ops.common import PREFERRED_PACKAGE_CLI_NAME

from .cli_commands import register_all_commands
from .cli_commands.common import context_from_args
from .container import run_in_docker


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=PREFERRED_PACKAGE_CLI_NAME)
    subparsers = parser.add_subparsers(dest="command", required=True)
    register_all_commands(subparsers)
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if (
            getattr(args, "engine", "host") == "docker"
            and os.getenv("FEELPP_PKG_IN_CONTAINER") != "1"
        ):
            context = context_from_args(args)
            try:
                run_in_docker(
                    context,
                    argv=argv,
                    image=getattr(args, "container_image", None),
                    state_root=getattr(args, "docker_state_root", None),
                    dry_run=getattr(args, "dry_run", False),
                )
            except subprocess.CalledProcessError as exc:
                return exc.returncode or 1
            return 0
        return int(args.func(args))
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
