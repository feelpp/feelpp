from __future__ import annotations

import argparse

from ....oci.commands import register_catalog_image_commands


def register_image_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    image_parser = subparsers.add_parser("image", help="Generate Spack-backed Docker image inputs")
    image_subparsers = image_parser.add_subparsers(dest="spack_image_command", required=True)
    register_catalog_image_commands(image_subparsers, backend_filter="spack")
