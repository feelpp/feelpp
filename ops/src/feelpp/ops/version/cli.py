from __future__ import annotations

import argparse
import json
import sys

from feelpp.ops.common import PREFERRED_VERSION_CLI_NAME
from feelpp.pkg.config import discover_repo_root

from .models import SemanticVersion
from .release import ReleaseService
from .repository import VersionRepository


def command_show(args: argparse.Namespace) -> int:
    repository = VersionRepository(repo_root=args.repo_root)
    print(json.dumps(repository.read_state().as_dict(), indent=2))
    return 0


def command_bump(args: argparse.Namespace) -> int:
    repository = VersionRepository(repo_root=args.repo_root)
    state = repository.bump_upstream(
        SemanticVersion.parse(args.version),
        message=args.message,
        dry_run=args.dry_run,
    )
    print(json.dumps(state.as_dict(), indent=2))
    return 0


def command_revision_bump(args: argparse.Namespace) -> int:
    repository = VersionRepository(repo_root=args.repo_root)
    state = repository.bump_revision(
        message=args.message,
        dists=tuple(args.dist) or None,
        dry_run=args.dry_run,
    )
    print(json.dumps(state.as_dict(), indent=2))
    return 0


def command_release(args: argparse.Namespace) -> int:
    service = ReleaseService(repo_root=args.repo_root)
    plan = service.execute_release(args.version, dry_run=args.dry_run, dists=tuple(args.dist) or None)
    print(json.dumps(plan.as_dict(), indent=2))
    return 0


def command_sync(args: argparse.Namespace) -> int:
    repository = VersionRepository(repo_root=args.repo_root)
    state = repository.sync_changelogs(
        message=args.message,
        dists=tuple(args.dist) or None,
        dry_run=args.dry_run,
    )
    print(json.dumps(state.as_dict(), indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=PREFERRED_VERSION_CLI_NAME)
    parser.add_argument(
        "--repo-root",
        default=str(discover_repo_root()),
        help="Path to the Feel++ repository root",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    show_parser = subparsers.add_parser("show", help="Show the current version state")
    show_parser.set_defaults(func=command_show)

    bump_parser = subparsers.add_parser("bump", help="Bump the upstream semantic version")
    bump_parser.add_argument("version", help="Semantic version to apply")
    bump_parser.add_argument(
        "--message",
        default="New upstream release",
        help="Debian changelog entry message",
    )
    bump_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview the updated version state without writing files",
    )
    bump_parser.set_defaults(func=command_bump)

    revision_parser = subparsers.add_parser("revision", help="Packaging revision commands")
    revision_subparsers = revision_parser.add_subparsers(dest="revision_command", required=True)
    revision_bump_parser = revision_subparsers.add_parser(
        "bump",
        help="Increment the Debian packaging revision without changing the upstream semantic version",
    )
    revision_bump_parser.add_argument(
        "--message",
        default="Packaging revision update",
        help="Debian changelog entry message",
    )
    revision_bump_parser.add_argument(
        "--dist",
        action="append",
        default=[],
        help="Restrict the revision bump to a distro. Repeat to target multiple distros.",
    )
    revision_bump_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview the updated version state without writing files",
    )
    revision_bump_parser.set_defaults(func=command_revision_bump)

    sync_parser = subparsers.add_parser(
        "sync",
        help="Sync Debian changelog heads from the root CMake version and manifest package revisions",
    )
    sync_parser.add_argument(
        "--message",
        default="Packaging metadata sync",
        help="Debian changelog entry message",
    )
    sync_parser.add_argument(
        "--dist",
        action="append",
        default=[],
        help="Restrict changelog sync to a distro. Repeat to target multiple distros.",
    )
    sync_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview the synced version state without writing files",
    )
    sync_parser.set_defaults(func=command_sync)

    release_parser = subparsers.add_parser("release", help="Create the git tag and GitHub release")
    release_parser.add_argument(
        "version",
        help="Upstream semantic version for GitHub releases or full Debian package version for packaging-only releases",
    )
    release_parser.add_argument(
        "--dist",
        action="append",
        default=[],
        help="Restrict release validation to a distro. Repeat to target multiple distros.",
    )
    release_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run all validations and preview the generated release content without publishing",
    )
    release_parser.set_defaults(func=command_release)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
