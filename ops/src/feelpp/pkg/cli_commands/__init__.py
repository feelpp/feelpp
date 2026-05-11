from __future__ import annotations

import argparse

from .build import register_build_commands
from .image import register_image_commands
from .inspect import register_inspect_commands
from .job import register_job_commands
from .pbuilder import register_pbuilder_commands
from .publish import register_publish_commands
from .repo import register_repo_commands


def register_all_commands(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    register_job_commands(subparsers)
    register_inspect_commands(subparsers)
    register_pbuilder_commands(subparsers)
    register_build_commands(subparsers)
    register_repo_commands(subparsers)
    register_publish_commands(subparsers)
    register_image_commands(subparsers)
