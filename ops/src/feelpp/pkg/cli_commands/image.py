from __future__ import annotations

import argparse

from ..backends.debian.commands.image import register_runtime_and_publish_commands
from ..oci.commands import register_catalog_image_commands


def register_image_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    image_parser = subparsers.add_parser("image", help="OCI image planning and runtime image commands")
    image_subparsers = image_parser.add_subparsers(dest="image_command", required=True)
    register_catalog_image_commands(image_subparsers)
    register_runtime_and_publish_commands(image_subparsers)
