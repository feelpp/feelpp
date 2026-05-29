from __future__ import annotations

import argparse

from ..backends.debian.commands import (
    register_backend_commands as register_debian_backend_commands,
    register_commands as register_debian_commands,
)
from ..backends.spack.commands import register_commands as register_spack_commands
from .image import register_image_commands as register_top_level_image_commands


BACKEND_COMMAND_REGISTRARS = (
    register_spack_commands,
    register_debian_backend_commands,
)


def register_all_commands(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    for registrar in BACKEND_COMMAND_REGISTRARS:
        registrar(subparsers)
    # Top-level Debian commands remain available as compatibility aliases.
    register_debian_commands(subparsers, include_image=False)
    register_top_level_image_commands(subparsers)
