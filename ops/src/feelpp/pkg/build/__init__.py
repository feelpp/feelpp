from __future__ import annotations

from pathlib import Path
import shutil

from ..config import PackagingContext
from ..debian import collect_component_internal_build_dependencies
from ..graph import BuildPlan
from ..localrepo import stage_outputs
from ..pbuilder import prepare_base, prepare_runtime_assets
from ..workspace import update_job_manifest_state, write_job_manifest
from . import archive as archive_module
from . import outer_prefix as outer_prefix_module
from . import sourcepkg as sourcepkg_module
from .elfcheck import validate_runtime_linkage
from .runner import run_pbuilder_build

__all__ = [
    "archive_module",
    "build_chain",
    "build_component",
    "clear_result_dir",
    "collect_component_internal_build_dependencies",
    "outer_prefix_module",
    "run_pbuilder_build",
    "sourcepkg_module",
    "stage_outputs",
    "validate_runtime_linkage",
]


def clear_result_dir(path: Path) -> None:
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True, exist_ok=True)


def build_component(
    context: PackagingContext,
    component: str,
    *,
    dry_run: bool = False,
    skip_pbuilder_prepare: bool = False,
) -> Path:
    prepare_runtime_assets(context)
    if not skip_pbuilder_prepare:
        prepare_base(context, dry_run=dry_run)

    result_dir = context.component_results_dir(component)
    clear_result_dir(result_dir)

    outer_prefix_module._bootstrap_outer_build_deps(
        context,
        component,
        dry_run=dry_run,
    )
    build_env = outer_prefix_module._outer_internal_env(context)
    archive_path = archive_module._build_source_archive(
        context,
        component,
        build_env=build_env,
        dry_run=dry_run,
    )
    _, dsc_path, _ = sourcepkg_module._prepare_source_tree(
        context,
        component,
        archive_path,
        dry_run=dry_run,
    )

    allow_public_fallback = not collect_component_internal_build_dependencies(
        context,
        component,
    )
    run_pbuilder_build(
        context,
        dsc_path,
        result_dir,
        allow_public_fallback=allow_public_fallback,
        dry_run=dry_run,
    )

    if not dry_run:
        validate_runtime_linkage(result_dir)
        stage_outputs(context, result_dir)
        update_job_manifest_state(
            context,
            state="built",
            built_component=component,
        )

    return result_dir


def build_chain(
    context: PackagingContext,
    plan: BuildPlan,
    *,
    dry_run: bool = False,
    prepare_pbuilder: bool = True,
) -> None:
    write_job_manifest(context, state="initialized", plan=plan)
    if prepare_pbuilder:
        prepare_base(context, dry_run=dry_run)

    for component in plan.components:
        build_component(
            context,
            component.name,
            dry_run=dry_run,
            skip_pbuilder_prepare=True,
        )
