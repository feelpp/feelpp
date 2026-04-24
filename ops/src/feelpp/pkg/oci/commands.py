from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from .apt import DEFAULT_VARIANT, describe_apt_target, generate_apt_bake
from .catalog import get_image_target, list_image_targets
from .common import add_workspace_arguments, workspace_from_args
from .spack import DEFAULT_SPACK_REF, describe_spack_target, generate_spack_bake
from ..shell import format_command, run_checked


BUILD_ARG_ENV_KEYS = (
    "FEELPP_GITHUB_TOKEN",
    "FEELPP_GIRDER_API_KEY",
    "FEELPP_CKAN_API_KEY",
    "FEELPP_CKAN_URL",
    "FEELPP_CKAN_ORGANIZATION",
)


def _resolved_backend_filter(args: argparse.Namespace) -> str | None:
    return getattr(args, "catalog_backend_filter", None) or getattr(args, "backend", None)


def _generate_image_metadata(args: argparse.Namespace) -> dict[str, object]:
    workspace = workspace_from_args(args)
    backend_filter = _resolved_backend_filter(args)
    target = get_image_target(
        args.target,
        repo_root=workspace.repo_root,
        profile=args.profile,
        backend=backend_filter,
    )

    if target.image_backend == "apt":
        if getattr(args, "tag", None):
            raise ValueError("--tag is only supported for Spack image generation")
        return generate_apt_bake(
            workspace=workspace,
            target=target,
            variant_name=getattr(args, "variant", DEFAULT_VARIANT),
            cc=getattr(args, "cc", "clang"),
            cxx=getattr(args, "cxx", "clang++"),
            cmake_flags=getattr(args, "cmake_flags", ""),
            base_image=getattr(args, "base_image", None),
            registry=getattr(args, "registry", None),
            namespace=getattr(args, "namespace", None),
            bake_target=getattr(args, "bake_target", None),
            component=getattr(args, "component", None),
            from_image=getattr(args, "from_image", None),
        )
    if target.image_backend == "spack":
        return generate_spack_bake(
            workspace=workspace,
            target=target,
            environment_name=getattr(args, "environment", None),
            base_image=getattr(args, "base_image", None),
            spack_ref=getattr(args, "spack_ref", DEFAULT_SPACK_REF),
            bake_target=getattr(args, "bake_target", None),
            image_tag=getattr(args, "tag", None),
            registry=getattr(args, "registry", None),
            namespace=getattr(args, "namespace", None),
            component=getattr(args, "component", None),
            from_image=getattr(args, "from_image", None),
            cc=getattr(args, "cc", "clang"),
            cxx=getattr(args, "cxx", "clang++"),
            cmake_flags=getattr(args, "cmake_flags", ""),
            spack_build_jobs=getattr(args, "spack_build_jobs", None),
            spack_concurrent_packages=getattr(args, "spack_concurrent_packages", None),
        )
    raise ValueError(f"Unsupported image backend: {target.image_backend}")


def _selected_bake_groups(
    args: argparse.Namespace,
    metadata: dict[str, object],
) -> list[str]:
    requested_groups = list(getattr(args, "group", []) or [])
    if requested_groups:
        return requested_groups
    groups = metadata.get("recommended_groups")
    if not isinstance(groups, list) or not groups or not all(isinstance(item, str) and item for item in groups):
        raise ValueError("image metadata must provide a non-empty recommended_groups list")
    return groups


def _build_arg_set_flags() -> tuple[list[str], list[str]]:
    actual: list[str] = []
    display: list[str] = []
    for key in BUILD_ARG_ENV_KEYS:
        value = os.getenv(key)
        if not value:
            continue
        actual.extend(["--set", f"*.args.{key}={value}"])
        display.extend(["--set", f"*.args.{key}=***"])
    return actual, display


def _add_image_generation_arguments(
    parser: argparse.ArgumentParser,
    *,
    backend_filter: str | None,
) -> None:
    add_workspace_arguments(parser)
    parser.add_argument(
        "--profile",
        default="images",
        help="Planning profile to inspect, defaults to images",
    )
    parser.add_argument(
        "--target",
        required=True,
        help="Image target from .github/plan-ci.json, for example ubuntu:noble or spack:openmpi",
    )
    parser.add_argument(
        "--base-image",
        default=None,
        help="Override the base image declared by the image target",
    )
    parser.add_argument(
        "--component",
        default=None,
        help="Generate bake metadata for a specific component path: env, feelpp, toolboxes, mor, or full",
    )
    parser.add_argument(
        "--from-image",
        default=None,
        help="Use an existing OCI image as the starting image for the selected component/full build",
    )
    parser.add_argument(
        "--bake-target",
        default=None,
        help="Override the generated docker buildx bake target root name",
    )
    parser.add_argument(
        "--registry",
        default=None,
        help="Override the OCI registry host, defaults to FEELPP_PKG_OCI_REGISTRY or ghcr.io",
    )
    parser.add_argument(
        "--namespace",
        default=None,
        help="Override the OCI namespace path, defaults to the namespace part of FEELPP_PKG_OCI_REPOSITORY",
    )
    if backend_filter is None or backend_filter == "apt":
        parser.add_argument(
            "--variant",
            default=DEFAULT_VARIANT,
            help="Docker metadata variant name for apt-backed targets",
        )
        parser.add_argument(
            "--cc",
            default="clang",
            help="C compiler used in generated apt-backed component builds",
        )
        parser.add_argument(
            "--cxx",
            default="clang++",
            help="C++ compiler used in generated apt-backed component builds",
        )
        parser.add_argument(
            "--cmake-flags",
            default="",
            help="Extra CMake flags passed into generated apt-backed component builds",
        )
    if backend_filter is None or backend_filter == "spack":
        parser.add_argument(
            "--environment",
            default=None,
            help="Override the repository-owned Spack environment name, for example cpu/openmpi",
        )
        parser.add_argument(
            "--spack-ref",
            default=DEFAULT_SPACK_REF,
            help="Git ref used when cloning Spack inside the generated image",
        )
        parser.add_argument(
            "--tag",
            default=None,
            help="Override the generated image tag for Spack-backed targets",
        )
        parser.add_argument(
            "--spack-build-jobs",
            type=int,
            default=None,
            help="Override Spack package-level build parallelism for generated Spack environment images",
        )
        parser.add_argument(
            "--spack-concurrent-packages",
            type=int,
            default=None,
            help="Override Spack install-level package concurrency for generated Spack environment images",
        )


def command_image_targets(args: argparse.Namespace) -> int:
    workspace = workspace_from_args(args)
    backend_filter = _resolved_backend_filter(args)
    targets = list_image_targets(workspace.repo_root, profile=args.profile, backend=backend_filter)
    described_targets = []
    for target in targets:
        if target.image_backend == "apt":
            described_targets.append(describe_apt_target(workspace.repo_root, target))
        elif target.image_backend == "spack":
            described_targets.append(describe_spack_target(workspace.repo_root, target))
        else:
            described_targets.append({**target.as_dict(), "supported": False})
    print(
        json.dumps(
            {
                "repo_root": str(workspace.repo_root),
                "profile": args.profile,
                "targets": described_targets,
            },
            indent=2,
        )
    )
    return 0


def command_image_bake(args: argparse.Namespace) -> int:
    metadata = _generate_image_metadata(args)
    print(json.dumps(metadata, indent=2))
    return 0


def command_image_build(args: argparse.Namespace) -> int:
    metadata = _generate_image_metadata(args)
    bake_groups = _selected_bake_groups(args, metadata)
    build_arg_flags, display_build_arg_flags = _build_arg_set_flags()

    command = ["docker", "buildx", "bake", "-f", str(metadata["bake_file"])]
    if args.push:
        command.append("--push")
    else:
        command.append("--load")
    command.extend(build_arg_flags)
    command.extend(bake_groups)

    if args.dry_run:
        display_command = ["docker", "buildx", "bake", "-f", str(metadata["bake_file"])]
        if args.push:
            display_command.append("--push")
        else:
            display_command.append("--load")
        display_command.extend(display_build_arg_flags)
        display_command.extend(bake_groups)
        print(format_command(display_command))
        return 0

    run_checked(command, cwd=Path(str(metadata["context_dir"])))
    return 0


def register_catalog_image_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    *,
    backend_filter: str | None = None,
) -> None:
    targets_parser = subparsers.add_parser(
        "targets",
        help="List OCI image targets from .github/plan-ci.json",
    )
    add_workspace_arguments(targets_parser)
    targets_parser.add_argument(
        "--profile",
        default="images",
        help="Planning profile to inspect, defaults to images",
    )
    if backend_filter is None:
        targets_parser.add_argument(
            "--backend",
            choices=["apt", "spack"],
            default=None,
            help="Restrict the output to a single image backend",
        )
    targets_parser.set_defaults(
        func=command_image_targets,
        catalog_backend_filter=backend_filter,
    )

    bake_parser = subparsers.add_parser(
        "bake",
        help="Generate a bake-ready OCI image context for a planned image target",
    )
    _add_image_generation_arguments(bake_parser, backend_filter=backend_filter)
    bake_parser.set_defaults(
        func=command_image_bake,
        catalog_backend_filter=backend_filter,
    )

    build_parser = subparsers.add_parser(
        "build",
        help="Generate and build OCI images for a planned image target",
    )
    _add_image_generation_arguments(build_parser, backend_filter=backend_filter)
    build_parser.add_argument(
        "--group",
        action="append",
        default=[],
        help="Override the generated recommended bake groups. Repeat to add more groups.",
    )
    build_parser.add_argument(
        "--push",
        action="store_true",
        help="Push built images to the configured OCI registry instead of loading them locally",
    )
    build_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the docker buildx bake command without executing it",
    )
    build_parser.set_defaults(
        func=command_image_build,
        catalog_backend_filter=backend_filter,
    )
