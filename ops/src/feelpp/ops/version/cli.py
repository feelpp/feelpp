from __future__ import annotations

import argparse
import json
import sys

from feelpp.ops.common import PREFERRED_VERSION_CLI_NAME
from feelpp.pkg.config import discover_repo_root

from .models import SemanticVersion
from .publications import HalPublicationService
from .release import ReleaseService
from .repository import VersionRepository


def _format_release_plan(plan) -> str:
    lines = [
        f"Release: {plan.title}",
        f"Tag: {plan.tag}",
        f"Kind: {plan.release_kind}",
        f"Branch: {plan.branch}",
        f"Channel: {plan.channel}",
        f"Head: {plan.head_sha}",
        f"Previous tag: {plan.previous_tag or '<none>'}",
        f"Prerelease: {'yes' if plan.prerelease else 'no'}",
        f"Dry run: {'yes' if plan.dry_run else 'no'}",
    ]

    lines.extend(["", "Package notes:", plan.package_notes])

    if plan.generated_notes_preview:
        lines.extend(["", "Generated notes preview:", plan.generated_notes_preview])

    return "\n".join(lines)


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
    plan = service.execute_release(
        args.version,
        dry_run=args.dry_run,
        dists=tuple(args.dist) or None,
        publication_rows=args.publications_rows,
        publication_since=args.publications_since,
    )
    if args.pretty:
        print(_format_release_plan(plan))
    else:
        print(json.dumps(plan.as_dict(), indent=2))
    return 0


def command_publications(args: argparse.Namespace) -> int:
    service = HalPublicationService(repo_root=args.repo_root)
    publications = service.fetch(
        query=args.query,
        rows=args.rows,
        sort=args.sort,
        collections=tuple(args.collection) or None,
        since=args.since,
    )
    if args.pretty:
        print(service.format_markdown(publications))
    else:
        print(json.dumps([publication.as_dict() for publication in publications], indent=2))
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
    release_parser.add_argument(
        "--pretty",
        action="store_true",
        help="Print the release plan in a human-readable format instead of JSON.",
    )
    release_parser.add_argument(
        "--publications-rows",
        type=int,
        default=None,
        help="Limit the number of HAL publications included in the release notes.",
    )
    release_parser.add_argument(
        "--publications-since",
        default=None,
        help="Limit HAL publications in the release notes to items produced since YYYY-MM-DD or a full RFC3339 timestamp.",
    )
    release_parser.set_defaults(func=command_release)

    publications_parser = subparsers.add_parser(
        "publications",
        help="Fetch recent Feel++ publications from HAL",
    )
    publications_parser.add_argument(
        "--query",
        default=None,
        help="Optional HAL query override. Defaults to the configured collection-driven query.",
    )
    publications_parser.add_argument(
        "--collection",
        action="append",
        default=[],
        help="Restrict HAL harvesting to a collection code. Repeat to target multiple collections.",
    )
    publications_parser.add_argument(
        "--rows",
        type=int,
        default=None,
        help="Maximum number of HAL records to return.",
    )
    publications_parser.add_argument(
        "--sort",
        default=None,
        help="HAL sort expression, for example 'producedDate_tdate desc'.",
    )
    publications_parser.add_argument(
        "--since",
        default=None,
        help="Optional lower bound for producedDate_tdate, as YYYY-MM-DD or full RFC3339 timestamp.",
    )
    publications_parser.add_argument(
        "--pretty",
        action="store_true",
        help="Render the publications as a Markdown release-notes section.",
    )
    publications_parser.set_defaults(func=command_publications)
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
