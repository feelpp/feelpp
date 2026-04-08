from __future__ import annotations

from pathlib import Path
import os
import re

from ..config import PackagingContext
from .common import (
    RepoPackageRef,
    aptly_base_command,
    build_publish_args,
    create_passphrase_file,
    publish_identity,
    publish_snapshot_from_repo,
    repo_name,
    run_capture,
    run_checked,
    run_probe,
    snapshot_id,
)


PRERELEASE_VERSION_RE = re.compile(
    r"^(?P<base>[^~]+)~(?P<stage>preview|alpha|beta|rc)(?:[.-].+)?$"
)


def parse_repo_search_output(output: str) -> list[RepoPackageRef]:
    refs: list[RepoPackageRef] = []
    for line in output.splitlines():
        raw = line.strip()
        if not raw:
            continue
        parts = raw.split("_", 2)
        if len(parts) != 3:
            continue
        refs.append(
            RepoPackageRef(
                raw=raw,
                name=parts[0],
                version=parts[1],
                architecture=parts[2],
            )
        )
    return refs


def _upstream_version(version: str) -> str:
    if "-" not in version:
        return version
    return version.rsplit("-", 1)[0]


def _prerelease_base(version: str) -> str | None:
    match = PRERELEASE_VERSION_RE.match(_upstream_version(version))
    if not match:
        return None
    return match.group("base")


def cleanup_candidates(package_refs: list[RepoPackageRef]) -> list[RepoPackageRef]:
    released_bases = {
        _upstream_version(package_ref.version)
        for package_ref in package_refs
        if _prerelease_base(package_ref.version) is None
    }
    return [
        package_ref
        for package_ref in package_refs
        if (base := _prerelease_base(package_ref.version)) is not None and base in released_bases
    ]


def _chunks(items: list[str], size: int) -> list[list[str]]:
    return [items[index : index + size] for index in range(0, len(items), size)]


def publish_cleanup(
    context: PackagingContext,
    *,
    dry_run: bool = False,
) -> None:
    env = context.shell_env()
    repo_name_value = repo_name(context)
    publish_distribution, publish_component, publish_target = publish_identity(context)
    aptly = aptly_base_command()

    if not run_probe(
        aptly + ["repo", "show", repo_name_value],
        cwd=context.repo_root,
        env=env,
    ):
        raise RuntimeError(f"Local aptly repository not found: {repo_name_value}")

    package_refs = parse_repo_search_output(
        run_capture(
            aptly + ["repo", "search", repo_name_value],
            cwd=context.repo_root,
            env=env,
        )
    )
    stale_refs = cleanup_candidates(package_refs)
    if not stale_refs:
        print(f"No prerelease packages eligible for cleanup in {repo_name_value}")
        return

    print(
        f"Cleaning {len(stale_refs)} prerelease package refs in {repo_name_value} "
        "whose final release already exists:"
    )
    for package_ref in stale_refs:
        print(f"  {package_ref.raw}")

    passphrase_file = None if dry_run else create_passphrase_file(context)
    try:
        publish_args = build_publish_args(
            passphrase_file=passphrase_file
            if passphrase_file is not None
            else (Path("/tmp/aptly-passphrase") if dry_run and os.getenv("GPG_PASSPHRASE") else None)
        )
        for refs_chunk in _chunks([package_ref.raw for package_ref in stale_refs], 100):
            run_checked(
                aptly + ["repo", "remove", repo_name_value, *refs_chunk],
                cwd=context.repo_root,
                env=env,
                dry_run=dry_run,
            )

        snapshot_name = f"{repo_name_value}-snapshot-cleanup-{snapshot_id()}"
        publish_snapshot_from_repo(
            context,
            env=env,
            snapshot_name=snapshot_name,
            publish_distribution=publish_distribution,
            publish_component=publish_component,
            publish_target=publish_target,
            publish_args=publish_args,
            dry_run=dry_run,
        )
        print(f"Published cleanup snapshot: {snapshot_name}")
    finally:
        if passphrase_file is not None and passphrase_file.exists():
            passphrase_file.unlink()
