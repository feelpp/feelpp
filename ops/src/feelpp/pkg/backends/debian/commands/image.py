from __future__ import annotations

import argparse

from ....image import build_runtime_image, publish_apptainer_image, publish_docker_image
from .common import add_context_arguments, context_from_args


def command_image_runtime(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    build_runtime_image(
        context,
        feelpp_version=args.feelpp_version,
        tag=args.tag,
        packages=args.package or None,
        base_image=args.base_image,
        apt_key_file=args.apt_key_file,
        gpg_key=args.gpg_key,
        repo_cache_token=args.repo_cache_token,
        no_cache=args.no_cache,
        dry_run=args.dry_run,
    )
    return 0


def command_image_publish_docker(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    publish_docker_image(
        context,
        source_ref=args.source_ref,
        target_ref=args.target_ref,
        registry=args.registry,
        repository=args.repository,
        tag=args.tag,
        dry_run=args.dry_run,
    )
    return 0


def command_image_publish_apptainer(args: argparse.Namespace) -> int:
    context = context_from_args(args)
    publish_apptainer_image(
        context,
        source_ref=args.source_ref,
        target_ref=args.target_ref,
        registry=args.registry,
        repository=args.repository,
        tag=args.tag,
        output_path=args.output,
        dry_run=args.dry_run,
    )
    return 0


def register_runtime_and_publish_commands(
    image_subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    image_runtime_parser = image_subparsers.add_parser(
        "runtime",
        help="Build a runtime image that installs pinned Feel++ packages from apt.feelpp.org",
    )
    add_context_arguments(image_runtime_parser)
    image_runtime_parser.add_argument(
        "--feelpp-version",
        required=True,
        help="Exact package version to install, for example 0.111.0~preview.13-1",
    )
    image_runtime_parser.add_argument(
        "--tag",
        default=None,
        help="Docker image tag to build locally",
    )
    image_runtime_parser.add_argument(
        "--base-image",
        default=None,
        help="Base image to use, defaults to <flavor>:<dist>",
    )
    image_runtime_parser.add_argument(
        "--apt-key-file",
        default=None,
        help="Path to a prebuilt Feel++ archive keyring file to copy into the image",
    )
    image_runtime_parser.add_argument(
        "--gpg-key",
        default=None,
        help="Public key fingerprint or key ID to export from the local GnuPG keyring into the image",
    )
    image_runtime_parser.add_argument(
        "--repo-cache-token",
        default=None,
        help="Override the Docker cache-busting token for the Feel++ apt install layer",
    )
    image_runtime_parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Force docker build to ignore all cached layers",
    )
    image_runtime_parser.add_argument(
        "--package",
        action="append",
        default=[],
        help="Package to install in the image. Repeat to add more packages.",
    )
    image_runtime_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the docker build command without executing it",
    )
    image_runtime_parser.set_defaults(func=command_image_runtime)

    image_publish_parser = image_subparsers.add_parser("publish", help="Publish image artifacts")
    image_publish_subparsers = image_publish_parser.add_subparsers(
        dest="image_publish_command",
        required=True,
    )

    image_publish_docker_parser = image_publish_subparsers.add_parser(
        "docker",
        help="Tag and push a local Docker image to an OCI registry",
    )
    add_context_arguments(image_publish_docker_parser)
    image_publish_docker_parser.add_argument(
        "--source-ref",
        required=True,
        help="Local Docker image reference",
    )
    image_publish_docker_parser.add_argument(
        "--target-ref",
        default=None,
        help="Fully-qualified OCI target reference",
    )
    image_publish_docker_parser.add_argument(
        "--registry",
        default=None,
        help="Registry host, defaults to ghcr.io",
    )
    image_publish_docker_parser.add_argument(
        "--repository",
        default=None,
        help="OCI repository path, defaults to feelpp/feelpp",
    )
    image_publish_docker_parser.add_argument(
        "--tag",
        default=None,
        help="Target image tag, defaults to the source tag",
    )
    image_publish_docker_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    image_publish_docker_parser.set_defaults(func=command_image_publish_docker)

    image_publish_apptainer_parser = image_publish_subparsers.add_parser(
        "apptainer",
        help="Convert a Docker image to a SIF artifact and push it via ORAS",
    )
    add_context_arguments(image_publish_apptainer_parser)
    image_publish_apptainer_parser.add_argument(
        "--source-ref",
        required=True,
        help="Local Docker image reference used as the Apptainer source",
    )
    image_publish_apptainer_parser.add_argument(
        "--target-ref",
        default=None,
        help="OCI target reference without or with the oras:// prefix",
    )
    image_publish_apptainer_parser.add_argument(
        "--registry",
        default=None,
        help="Registry host, defaults to ghcr.io",
    )
    image_publish_apptainer_parser.add_argument(
        "--repository",
        default=None,
        help="OCI repository path, defaults to feelpp/feelpp",
    )
    image_publish_apptainer_parser.add_argument(
        "--tag",
        default=None,
        help="Target artifact tag, defaults to <source-tag>_sif",
    )
    image_publish_apptainer_parser.add_argument(
        "--output",
        default=None,
        help="Local path for the generated .sif file",
    )
    image_publish_apptainer_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    image_publish_apptainer_parser.set_defaults(func=command_image_publish_apptainer)


def register_image_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    image_parser = subparsers.add_parser("image", help="Docker image commands")
    image_subparsers = image_parser.add_subparsers(dest="image_command", required=True)
    register_runtime_and_publish_commands(image_subparsers)
