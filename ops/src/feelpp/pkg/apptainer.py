from __future__ import annotations

from pathlib import Path
import json
import os
import shutil

from .config import PackagingContext
from .docker import default_publish_ref, image_ref_tag, safe_workdir_name
from .shell import run
from .workspace import ensure_workspace


DEFAULT_APPTAINER_TAG_SUFFIX = "_sif"
DEFAULT_APPTAINER_ORAS_PREFIX = "oras://"
DEFAULT_APPTAINER_DOCKER_PREFIX = "docker://"
DEFAULT_APPTAINER_DOCKER_DAEMON_PREFIX = "docker-daemon:"
DEFAULT_APPTAINER_BINARY_ENV = "FEELPP_PKG_APPTAINER_BIN"
DEFAULT_APPTAINER_BINARY_FALLBACK = "/opt/apptainer/latest/bin/apptainer"


def default_apptainer_tag(source_ref: str, *, tag: str | None = None) -> str:
    if tag:
        return tag
    return f"{image_ref_tag(source_ref)}{DEFAULT_APPTAINER_TAG_SUFFIX}"


def default_apptainer_publish_ref(
    source_ref: str,
    *,
    target_ref: str | None = None,
    registry: str | None = None,
    repository: str | None = None,
    tag: str | None = None,
) -> str:
    return default_publish_ref(
        source_ref,
        target_ref=target_ref,
        registry=registry,
        repository=repository,
        tag=default_apptainer_tag(source_ref, tag=tag),
    )


def as_oras_ref(target_ref: str) -> str:
    if target_ref.startswith(DEFAULT_APPTAINER_ORAS_PREFIX):
        return target_ref
    return f"{DEFAULT_APPTAINER_ORAS_PREFIX}{target_ref}"


def as_apptainer_source_ref(source_ref: str) -> str:
    if source_ref.startswith(
        (
            DEFAULT_APPTAINER_DOCKER_PREFIX,
            DEFAULT_APPTAINER_DOCKER_DAEMON_PREFIX,
            DEFAULT_APPTAINER_ORAS_PREFIX,
            "library://",
        )
    ):
        return source_ref
    if "/" in source_ref:
        return f"{DEFAULT_APPTAINER_DOCKER_PREFIX}{source_ref}"
    return f"{DEFAULT_APPTAINER_DOCKER_DAEMON_PREFIX}{source_ref}"


def default_apptainer_binary() -> str:
    configured = os.getenv(DEFAULT_APPTAINER_BINARY_ENV)
    if configured:
        return configured
    resolved = shutil.which("apptainer")
    if resolved:
        return resolved
    if Path(DEFAULT_APPTAINER_BINARY_FALLBACK).is_file():
        return DEFAULT_APPTAINER_BINARY_FALLBACK
    versioned_candidates = sorted(
        Path("/opt/apptainer").glob("v*/apptainer/bin/apptainer"),
        reverse=True,
    )
    for candidate in versioned_candidates:
        if candidate.is_file():
            return str(candidate)
    return "apptainer"


def default_sif_path(
    context: PackagingContext,
    *,
    target_ref: str,
) -> Path:
    safe_target = safe_workdir_name(target_ref)
    output_dir = context.job_root / "images" / safe_target
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / f"{safe_target}.sif"


def publish_apptainer_image(
    context: PackagingContext,
    *,
    source_ref: str,
    target_ref: str | None = None,
    registry: str | None = None,
    repository: str | None = None,
    tag: str | None = None,
    output_path: str | None = None,
    dry_run: bool = False,
) -> dict[str, str]:
    ensure_workspace(context)
    resolved_target_ref = default_apptainer_publish_ref(
        source_ref,
        target_ref=target_ref,
        registry=registry,
        repository=repository,
        tag=tag,
    )
    resolved_oras_ref = as_oras_ref(resolved_target_ref)
    resolved_source_ref = as_apptainer_source_ref(source_ref)
    resolved_output_path = (
        Path(output_path).expanduser().resolve()
        if output_path
        else default_sif_path(context, target_ref=resolved_target_ref)
    )
    resolved_output_path.parent.mkdir(parents=True, exist_ok=True)
    apptainer_binary = default_apptainer_binary()

    run(
        [apptainer_binary, "build", str(resolved_output_path), resolved_source_ref],
        cwd=context.repo_root,
        env=context.shell_env(),
        dry_run=dry_run,
    )
    run(
        [apptainer_binary, "push", str(resolved_output_path), resolved_oras_ref],
        cwd=context.repo_root,
        env=context.shell_env(),
        dry_run=dry_run,
    )
    result = {
        "source_ref": source_ref,
        "resolved_source_ref": resolved_source_ref,
        "target_ref": resolved_target_ref,
        "oras_ref": resolved_oras_ref,
        "output_path": str(resolved_output_path),
        "apptainer_binary": apptainer_binary,
    }
    print(json.dumps(result, indent=2))
    return result
