from __future__ import annotations

import argparse
import json

from .apt import DEFAULT_VARIANT, describe_apt_target, generate_apt_bake
from .catalog import get_image_target, list_image_targets
from .common import add_workspace_arguments, workspace_from_args
from .spack import DEFAULT_SPACK_REF, describe_spack_target, generate_spack_bake


def _resolved_backend_filter(args: argparse.Namespace) -> str | None:
    return getattr(args, "catalog_backend_filter", None) or getattr(args, "backend", None)


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
        metadata = generate_apt_bake(
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
        )
    elif target.image_backend == "spack":
        metadata = generate_spack_bake(
            workspace=workspace,
            target=target,
            environment_name=getattr(args, "environment", None),
            base_image=getattr(args, "base_image", None),
            spack_ref=getattr(args, "spack_ref", DEFAULT_SPACK_REF),
            bake_target=getattr(args, "bake_target", None),
            image_tag=getattr(args, "tag", None),
            registry=getattr(args, "registry", None),
            namespace=getattr(args, "namespace", None),
        )
    else:
        raise ValueError(f"Unsupported image backend: {target.image_backend}")

    print(json.dumps(metadata, indent=2))
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
    add_workspace_arguments(bake_parser)
    bake_parser.add_argument(
        "--profile",
        default="images",
        help="Planning profile to inspect, defaults to images",
    )
    bake_parser.add_argument(
        "--target",
        required=True,
        help="Image target from .github/plan-ci.json, for example ubuntu:noble or spack:openmpi",
    )
    bake_parser.add_argument(
        "--base-image",
        default=None,
        help="Override the base image declared by the image target",
    )
    bake_parser.add_argument(
        "--bake-target",
        default=None,
        help="Override the generated docker buildx bake target root name",
    )
    bake_parser.add_argument(
        "--registry",
        default=None,
        help="Override the OCI registry host, defaults to FEELPP_PKG_OCI_REGISTRY or ghcr.io",
    )
    bake_parser.add_argument(
        "--namespace",
        default=None,
        help="Override the OCI namespace path, defaults to the namespace part of FEELPP_PKG_OCI_REPOSITORY",
    )
    if backend_filter is None or backend_filter == "apt":
        bake_parser.add_argument(
            "--variant",
            default=DEFAULT_VARIANT,
            help="Docker metadata variant name for apt-backed targets",
        )
        bake_parser.add_argument(
            "--cc",
            default="clang",
            help="C compiler used in generated apt-backed component builds",
        )
        bake_parser.add_argument(
            "--cxx",
            default="clang++",
            help="C++ compiler used in generated apt-backed component builds",
        )
        bake_parser.add_argument(
            "--cmake-flags",
            default="",
            help="Extra CMake flags passed into generated apt-backed component builds",
        )
    if backend_filter is None or backend_filter == "spack":
        bake_parser.add_argument(
            "--environment",
            default=None,
            help="Override the repository-owned Spack environment name, for example cpu/openmpi",
        )
        bake_parser.add_argument(
            "--spack-ref",
            default=DEFAULT_SPACK_REF,
            help="Git ref used when cloning Spack inside the generated image",
        )
        bake_parser.add_argument(
            "--tag",
            default=None,
            help="Override the generated image tag for Spack-backed targets",
        )
    bake_parser.set_defaults(
        func=command_image_bake,
        catalog_backend_filter=backend_filter,
    )
