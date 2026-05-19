from __future__ import annotations

import argparse

from .build import register_build_commands
from .image import register_image_commands
from .inspect import register_inspect_commands
from .job import register_job_commands
from .pbuilder import register_pbuilder_commands
from .publish import register_publish_commands
from .repo import register_repo_commands


DEBIAN_COMMAND_REGISTRARS = (
    register_job_commands,
    register_inspect_commands,
    register_pbuilder_commands,
    register_build_commands,
    register_repo_commands,
    register_publish_commands,
    register_image_commands,
)


def register_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    *,
    include_image: bool = True,
) -> None:
    for registrar in DEBIAN_COMMAND_REGISTRARS:
        if not include_image and registrar is register_image_commands:
            continue
        registrar(subparsers)


def register_backend_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    debian_parser = subparsers.add_parser("debian", help="Debian packaging commands")
    debian_subparsers = debian_parser.add_subparsers(dest="debian_area", required=True)
    register_commands(debian_subparsers, include_image=True)
