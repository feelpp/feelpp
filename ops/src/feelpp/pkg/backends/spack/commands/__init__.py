from __future__ import annotations

import argparse

from .env import register_env_commands
from .image import register_image_commands
from .local import register_local_commands


def register_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    spack_parser = subparsers.add_parser("spack", help="Spack commands")
    spack_subparsers = spack_parser.add_subparsers(dest="spack_command", required=True)
    register_local_commands(spack_subparsers)
    register_env_commands(spack_subparsers)
    register_image_commands(spack_subparsers)
