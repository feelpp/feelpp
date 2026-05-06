from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys

from feelpp.ops.common import PACKAGE_CLI_ALIAS_COMMANDS, PREFERRED_PACKAGE_CLI_NAME

from .cli_commands import register_all_commands
from .cli_commands.common import context_from_args
from .container import run_in_docker


def build_parser(*, prog: str | None = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=prog or PREFERRED_PACKAGE_CLI_NAME)
    subparsers = parser.add_subparsers(dest="command", required=True)
    register_all_commands(subparsers)
    return parser


def _normalize_alias_argv(invoked_name: str, argv: list[str]) -> list[str]:
    alias_command = PACKAGE_CLI_ALIAS_COMMANDS.get(invoked_name)
    if not alias_command:
        return argv
    if tuple(argv[: len(alias_command)]) == alias_command:
        return argv
    return [*alias_command, *argv]


def main(argv: list[str] | None = None) -> int:
    invoked_name = Path(sys.argv[0]).name
    argv = list(sys.argv[1:] if argv is None else argv)
    argv = _normalize_alias_argv(invoked_name, argv)
    parser_prog = invoked_name if invoked_name in PACKAGE_CLI_ALIAS_COMMANDS else None
    parser = build_parser(prog=parser_prog)
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
