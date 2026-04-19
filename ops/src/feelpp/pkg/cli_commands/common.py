from __future__ import annotations

import argparse

from ..config import PackagingContext
from ..graph import build_plan, load_manifest


def add_context_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--repo-root", help="Path to the Feel++ repository root")
    parser.add_argument("--dist", default=None, help="Distribution name, for example noble")
    parser.add_argument("--flavor", default=None, help="Flavor name, for example ubuntu")
    parser.add_argument("--branch", default=None, help="Branch name used to derive the channel")
    parser.add_argument("--channel", default=None, help="Package channel, for example latest or stable")
    parser.add_argument("--job-id", default=None, help="Packaging job identifier")
    parser.add_argument("--job-root", default=None, help="Packaging job root directory")


def add_runtime_arguments(
    parser: argparse.ArgumentParser,
    *,
    default_engine: str = "docker",
) -> None:
    parser.add_argument(
        "--engine",
        choices=["host", "docker"],
        default=default_engine,
        help="Execution backend",
    )
    parser.add_argument(
        "--container-image",
        default=None,
        help="Override the pkg-env image used with --engine docker",
    )
    parser.add_argument(
        "--docker-state-root",
        default=None,
        help="Host directory used for persistent docker-side chroots and caches",
    )


def context_from_args(args: argparse.Namespace) -> PackagingContext:
    return PackagingContext.create(
        repo_root=args.repo_root,
        dist=args.dist,
        flavor=args.flavor,
        branch=args.branch,
        channel=args.channel,
        job_id=args.job_id,
        job_root=args.job_root,
    )


def parse_component_list(raw: str | None) -> list[str] | None:
    if not raw:
        return None
    return [item.strip() for item in raw.split(",") if item.strip()]


def load_plan(args: argparse.Namespace, context: PackagingContext):
    manifest = load_manifest(context.manifest_path)
    return build_plan(
        manifest,
        dist=context.dist,
        requested_components=parse_component_list(args.components),
        skipped_components=set(args.skip_component or ()),
        skip_text=args.skip_text,
        publish=getattr(args, "publish", False),
    )
