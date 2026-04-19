from __future__ import annotations

from pathlib import Path

from ..config import PackagingContext
from ..shell import run


def _component_preset(component: str) -> str:
    return component[len("feelpp-") :] if component.startswith("feelpp-") else component


def _build_source_archive(
    context: PackagingContext,
    component: str,
    *,
    build_env: dict[str, str] | None = None,
    dry_run: bool = False,
) -> Path:
    preset = _component_preset(component)
    run(
        [
            "cmake",
            "--preset",
            preset,
            "-DFEELPP_ENABLE_GIT=OFF",
            "-DLIBBSON_DIR=/usr",
            "-DLIBMONGOC_DIR=/usr",
        ],
        cwd=context.repo_root,
        env=build_env,
        dry_run=dry_run,
    )
    run(
        ["cmake", "--build", "--preset", preset, "-t", "dist"],
        cwd=context.repo_root,
        env=build_env,
        dry_run=dry_run,
    )

    build_dir = context.repo_root / "build" / preset
    candidates = sorted(build_dir.glob(f"{component}-*.tar.gz"))
    if candidates:
        return candidates[-1]
    if dry_run:
        return build_dir / f"{component}-source.tar.gz"
    raise FileNotFoundError(f"Source archive not found for {component} in {build_dir}")
