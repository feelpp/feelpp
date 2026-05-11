from __future__ import annotations

from pathlib import Path
import json

from .config import PackagingContext
from .graph import BuildPlan


def ensure_workspace(context: PackagingContext) -> None:
    for path in (
        context.job_root,
        context.local_repo_dir,
        context.artifacts_dir,
        context.results_dir,
        context.pbuilder_runtime_hookdir,
        context.pbuilder_keyrings_dir,
    ):
        path.mkdir(parents=True, exist_ok=True)


def write_job_manifest(
    context: PackagingContext,
    *,
    state: str,
    plan: BuildPlan | None = None,
    built_components: list[str] | None = None,
) -> Path:
    ensure_workspace(context)
    payload = {
        "state": state,
        "context": context.as_dict(),
        "plan": plan.as_dict() if plan else None,
        "built_components": built_components or [],
    }
    context.job_manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return context.job_manifest_path


def read_job_manifest(context: PackagingContext) -> dict[str, object] | None:
    path = context.job_manifest_path
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def update_job_manifest_state(
    context: PackagingContext,
    *,
    state: str | None = None,
    built_component: str | None = None,
) -> Path:
    ensure_workspace(context)
    payload = read_job_manifest(context) or {
        "state": "initialized",
        "context": context.as_dict(),
        "plan": None,
        "built_components": [],
    }
    payload["context"] = context.as_dict()
    if state is not None:
        payload["state"] = state

    built_components = list(payload.get("built_components") or [])
    if built_component and built_component not in built_components:
        built_components.append(built_component)
    payload["built_components"] = built_components

    context.job_manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return context.job_manifest_path
