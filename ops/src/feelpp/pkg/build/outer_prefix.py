from __future__ import annotations

from pathlib import Path
import os
import subprocess

from ..config import PackagingContext
from ..debian import (
    collect_component_build_dependencies,
    collect_component_internal_build_dependencies,
)
from ..shell import run


def _outer_root(context: PackagingContext) -> Path:
    return context.job_root / "outer-prefix"


def _outer_usr(context: PackagingContext) -> Path:
    return _outer_root(context) / "usr"


def _outer_repo_list_path(context: PackagingContext) -> Path:
    return Path(
        os.getenv(
            "FEELPP_PKG_OUTER_LOCAL_REPO_LIST",
            str(_outer_root(context) / "etc" / "apt" / "sources.list.d" / "feelpp-local-repo.list"),
        )
    )


def _outer_repo_prefs_path(context: PackagingContext) -> Path:
    return Path(
        os.getenv(
            "FEELPP_PKG_OUTER_LOCAL_REPO_PREFS",
            str(_outer_root(context) / "etc" / "apt" / "preferences.d" / "feelpp-local-repo.pref"),
        )
    )


def _outer_apt_archives_dir(context: PackagingContext) -> Path:
    return Path(
        os.getenv(
            "FEELPP_PKG_OUTER_APT_ARCHIVES_DIR",
            str(_outer_root(context) / "var" / "cache" / "apt" / "archives"),
        )
    )


def _outer_internal_env(context: PackagingContext) -> dict[str, str]:
    prefix = _outer_usr(context)
    lib_dir = prefix / "lib" / "x86_64-linux-gnu"
    env = dict(os.environ)
    env["FEELPP_DIR"] = str(prefix)
    env["PATH"] = _prepend_env_path(env.get("PATH"), prefix / "bin")
    env["CMAKE_PREFIX_PATH"] = _prepend_env_path(env.get("CMAKE_PREFIX_PATH"), prefix)
    env["LD_LIBRARY_PATH"] = _prepend_env_path(env.get("LD_LIBRARY_PATH"), lib_dir)
    env["PKG_CONFIG_PATH"] = _prepend_multi_path(
        env.get("PKG_CONFIG_PATH"),
        [
            lib_dir / "pkgconfig",
            prefix / "share" / "pkgconfig",
        ],
    )
    env["PYTHONPATH"] = _prepend_env_path(
        env.get("PYTHONPATH"),
        prefix / "lib" / "python3" / "dist-packages",
    )
    return env


def _prepend_env_path(current: str | None, path: Path) -> str:
    value = str(path)
    if not current:
        return value
    return f"{value}:{current}"


def _prepend_multi_path(current: str | None, paths: list[Path]) -> str:
    values = [str(path) for path in paths]
    if current:
        values.append(current)
    return ":".join(values)


def _internal_repo_archive_candidates(context: PackagingContext) -> list[Path]:
    patterns = ("feelpp*.deb", "libfeelpp*.deb", "python3-feelpp*.deb")
    candidates: list[Path] = []
    for pattern in patterns:
        candidates.extend(sorted(context.local_repo_dir.glob(pattern)))
    unique: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        unique.append(candidate)
    return unique


def _write_local_repo_files(context: PackagingContext, internal_packages: list[str]) -> None:
    repo_list_path = _outer_repo_list_path(context)
    repo_list_path.parent.mkdir(parents=True, exist_ok=True)
    repo_list_path.write_text(
        f"deb [trusted=yes] file://{context.local_repo_dir} ./\n",
        encoding="utf-8",
    )

    prefs_path = _outer_repo_prefs_path(context)
    prefs_path.parent.mkdir(parents=True, exist_ok=True)
    packages = " ".join(internal_packages) if internal_packages else "feelpp-* libfeelpp* python3-feelpp*"
    prefs_path.write_text(
        "\n".join(
            [
                f"Package: {packages}",
                "Pin: origin apt.feelpp.org",
                "Pin-Priority: 100",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _clear_apt_archives(context: PackagingContext) -> None:
    archives_dir = _outer_apt_archives_dir(context)
    archives_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("feelpp*.deb", "libfeelpp*.deb", "python3-feelpp*.deb"):
        for archive in archives_dir.glob(pattern):
            archive.unlink()


def _bootstrap_outer_build_deps(
    context: PackagingContext,
    component: str,
    *,
    dry_run: bool = False,
) -> None:
    outer_root = _outer_root(context)
    outer_root.mkdir(parents=True, exist_ok=True)

    build_dependencies = collect_component_build_dependencies(context, component)
    internal_packages = collect_component_internal_build_dependencies(context, component)

    _write_local_repo_files(context, internal_packages)
    _clear_apt_archives(context)

    if internal_packages:
        for archive in _internal_repo_archive_candidates(context):
            subprocess.run(
                ["dpkg-deb", "-x", str(archive), str(outer_root)],
                check=True,
            )

    env = _outer_internal_env(context)
    run(["apt-get", "update"], env=env, dry_run=dry_run)
    if build_dependencies:
        run(
            ["apt-get", "install", "-y", *build_dependencies],
            env=env,
            dry_run=dry_run,
        )
