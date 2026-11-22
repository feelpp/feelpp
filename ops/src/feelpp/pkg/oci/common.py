# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re

from ..core.context import WorkspaceContext, discover_repo_root
from ..docker import default_oci_registry, default_oci_repository


def resolve_repo_root(raw: str | None) -> Path:
    return Path(raw).expanduser().resolve() if raw else discover_repo_root()


def workspace_from_args(args: argparse.Namespace) -> WorkspaceContext:
    return WorkspaceContext.create(
        repo_root=args.repo_root,
        branch=getattr(args, "branch", None),
        channel=getattr(args, "channel", None),
        job_id=getattr(args, "job_id", None),
        job_root=getattr(args, "job_root", None),
    )


def add_workspace_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--repo-root", help="Path to the Feel++ repository root")
    parser.add_argument("--branch", default=None, help="Branch name used for image metadata")
    parser.add_argument("--channel", default=None, help="Channel name used for image metadata")
    parser.add_argument("--job-id", default=None, help="Packaging job identifier")
    parser.add_argument("--job-root", default=None, help="Packaging job root directory")


def safe_name(raw: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_.-]+", "-", raw).strip("-")


def default_bake_target_name(target: str) -> str:
    return safe_name(target.replace(":", "-"))


def branch_tag_suffix(branch: str) -> str:
    normalized = safe_name(branch)
    if not normalized or normalized == "develop":
        return ""
    return f"-{normalized}"


def branch_image_tag(oci_dist: str, *, branch: str, suffix: str = "") -> str:
    """Keep existing tags, shortening only oversized branch names to fit OCI tags."""
    branch_suffix = branch_tag_suffix(branch)
    tag = f"{oci_dist}{branch_suffix}{suffix}"
    if len(tag) <= 128:
        return tag

    # Retain the distribution and image role, and distinguish branches whose
    # readable prefixes are identical. Hash the original, unsanitized name.
    digest = hashlib.sha256(branch.encode("utf-8")).hexdigest()[:12]
    branch_budget = 128 - len(oci_dist) - len(suffix) - len(digest) - 2
    if not branch_suffix or branch_budget < 1:
        raise ValueError("OCI tag distribution and suffix leave no room for a branch name")
    shortened = branch_suffix[1:][:branch_budget].rstrip("-.")
    return f"{oci_dist}-{shortened}-{digest}{suffix}"


def default_oci_namespace() -> str:
    repository = default_oci_repository()
    parts = [part for part in repository.split("/") if part]
    if len(parts) > 1:
        return "/".join(parts[:-1])
    return parts[0] if parts else "feelpp"


def oci_image_ref(
    image_name: str,
    tag: str,
    *,
    registry: str | None = None,
    namespace: str | None = None,
) -> str:
    resolved_registry = registry or default_oci_registry()
    resolved_namespace = namespace or default_oci_namespace()
    repository = f"{resolved_namespace}/{image_name}" if resolved_namespace else image_name
    return f"{resolved_registry}/{repository}:{tag}"
