from __future__ import annotations

from pathlib import Path
import os

from ..config import PackagingContext
from .constants import CONTAINER_EXTRA_ROOT, CONTAINER_JOB_ROOT, CONTAINER_REPO_ROOT


def _strip_host_only_args(argv: list[str]) -> list[str]:
    stripped: list[str] = []
    skip_next = False
    host_only_flags = {
        "--engine",
        "--container-image",
        "--docker-state-root",
        "--repo-root",
        "--job-root",
    }
    for token in argv:
        if skip_next:
            skip_next = False
            continue
        if token in host_only_flags:
            skip_next = True
            continue
        if token.startswith("--engine="):
            continue
        if token.startswith("--container-image="):
            continue
        if token.startswith("--docker-state-root="):
            continue
        if token.startswith("--repo-root="):
            continue
        if token.startswith("--job-root="):
            continue
        stripped.append(token)
    return stripped


def _host_path_from_arg(raw: str, *, repo_root: Path) -> Path:
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = (repo_root / path).resolve()
    else:
        path = path.resolve()
    return path


def _container_path_for_host_path(
    host_path: Path,
    *,
    context: PackagingContext,
    extra_mounts: list[tuple[Path, Path]],
    fallback_name: str,
) -> Path:
    if host_path == context.repo_root or context.repo_root in host_path.parents:
        return CONTAINER_REPO_ROOT / host_path.relative_to(context.repo_root)
    if host_path == context.job_root or context.job_root in host_path.parents:
        return CONTAINER_JOB_ROOT / host_path.relative_to(context.job_root)

    container_path = CONTAINER_EXTRA_ROOT / fallback_name
    extra_mounts.append((host_path, container_path))
    return container_path


def _containerized_argv(context: PackagingContext, argv: list[str]) -> tuple[list[str], list[tuple[Path, Path]]]:
    inner = _strip_host_only_args(argv)
    extra_mounts: list[tuple[Path, Path]] = []

    index = 0
    while index < len(inner):
        token = inner[index]
        if token == "--input-dir" and index + 1 < len(inner):
            host_path = _host_path_from_arg(inner[index + 1], repo_root=context.repo_root)
            container_path = _container_path_for_host_path(
                host_path,
                context=context,
                extra_mounts=extra_mounts,
                fallback_name="publish-input",
            )
            inner[index + 1] = str(container_path)
            index += 2
            continue
        if token.startswith("--input-dir="):
            host_path = _host_path_from_arg(token.split("=", 1)[1], repo_root=context.repo_root)
            container_path = _container_path_for_host_path(
                host_path,
                context=context,
                extra_mounts=extra_mounts,
                fallback_name="publish-input",
            )
            inner[index] = f"--input-dir={container_path}"
        index += 1

    if len(inner) >= 3 and inner[0] == "repo" and inner[1] == "stage":
        host_path = _host_path_from_arg(inner[2], repo_root=context.repo_root)
        container_path = _container_path_for_host_path(
            host_path,
            context=context,
            extra_mounts=extra_mounts,
            fallback_name="result-dir",
        )
        inner[2] = str(container_path)

    return [
        *inner,
        "--repo-root",
        str(CONTAINER_REPO_ROOT),
        "--job-root",
        str(CONTAINER_JOB_ROOT),
    ], extra_mounts


def publish_uses_signing(argv: list[str]) -> bool:
    if os.getenv("FEELPP_APTLY_SKIP_SIGNING", "").lower() == "true":
        return False

    inner = _strip_host_only_args(argv)
    if len(inner) >= 2 and inner[0] == "publish" and inner[1] in {"snapshot", "switch", "cleanup"}:
        return True
    return len(inner) >= 2 and inner[0] == "build" and inner[1] == "chain" and "--publish" in inner


def publish_uses_aptly(argv: list[str]) -> bool:
    inner = _strip_host_only_args(argv)
    if len(inner) >= 2 and inner[0] == "publish" and inner[1] in {"snapshot", "switch", "cleanup"}:
        return True
    return len(inner) >= 2 and inner[0] == "build" and inner[1] == "chain" and "--publish" in inner


def _pass_through_env() -> set[str]:
    allowed = {
        "DEB_BUILD_OPTIONS",
        "GITHUB_RUN_ID",
        "GITHUB_RUN_ATTEMPT",
        "GITHUB_JOB",
        "GITHUB_REF_NAME",
        "BUILDKITE_BRANCH",
        "BUILDKITE_AGENT_NAME",
        "AWS_ACCESS_KEY_ID",
        "AWS_SECRET_ACCESS_KEY",
        "AWS_SESSION_TOKEN",
        "GPG_KEY",
        "GPG_PASSPHRASE",
        "FEELPP_APTLY_PUBLISH_TARGET",
        "FEELPP_APTLY_PUBLISH_DISTRIBUTION",
        "FEELPP_APTLY_PUBLISH_COMPONENT",
        "FEELPP_APTLY_SKIP_SIGNING",
        "FEELPP_PKG_SNAPSHOT_ID",
    }
    return {key for key in allowed if os.getenv(key)}
