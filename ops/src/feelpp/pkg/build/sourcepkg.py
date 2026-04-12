from __future__ import annotations

from pathlib import Path
import shutil
import tarfile

from ..config import PackagingContext
from ..shell import run


def _archive_version(component: str, archive_path: Path) -> str:
    stem = archive_path.name
    if stem.endswith(".tar.gz"):
        stem = stem[: -len(".tar.gz")]
    prefix = f"{component}-"
    if not stem.startswith(prefix):
        raise ValueError(f"Unexpected source archive name: {archive_path.name}")
    raw_version = stem[len(prefix) :]
    head, separator, tail = raw_version.partition("-")
    return head if not separator else f"{head}~{tail}"


def _orig_archive_name(component: str, version: str) -> str:
    return f"{component}_{version}.orig.tar.gz"


def _archive_root_members(handle: tarfile.TarFile) -> set[str]:
    roots: set[str] = set()
    for member in handle.getmembers():
        parts = Path(member.name).parts
        if not parts:
            continue
        if parts[0] == ".":
            if len(parts) == 1:
                continue
            roots.add(parts[1])
        else:
            roots.add(parts[0])
    return roots


def _prepare_source_tree(
    context: PackagingContext,
    component: str,
    archive_path: Path,
    *,
    dry_run: bool = False,
) -> tuple[Path, Path, str]:
    version = _archive_version(component, archive_path)
    source_root = context.job_root / "source-packages" / component
    tree_root = source_root / f"{component}-{version}"
    packaging_dir = (
        context.repo_root / "packaging" / "debian" / component / context.dist / "debian"
    )
    dsc_path = source_root / f"{component}_{version}-1.dsc"

    if dry_run:
        return tree_root, dsc_path, version

    shutil.rmtree(source_root, ignore_errors=True)
    source_root.mkdir(parents=True, exist_ok=True)

    orig_archive = source_root / _orig_archive_name(component, version)
    shutil.copy2(archive_path, orig_archive)
    with tarfile.open(orig_archive, "r:gz") as handle:
        archive_roots = _archive_root_members(handle)
        handle.extractall(source_root)

    if len(archive_roots) != 1:
        roots_text = ", ".join(sorted(archive_roots)) or "<none>"
        raise ValueError(
            f"Expected a single top-level directory in {orig_archive.name}, found {roots_text}"
        )

    extracted_root = source_root / next(iter(archive_roots))
    if not extracted_root.exists():
        raise FileNotFoundError(f"Extracted source root not found: {extracted_root}")
    if extracted_root != tree_root:
        extracted_root.rename(tree_root)

    if not packaging_dir.is_dir():
        raise FileNotFoundError(f"Packaging tree not found: {packaging_dir}")
    shutil.copytree(packaging_dir, tree_root / "debian", dirs_exist_ok=True)

    run(
        [
            "dch",
            "-v",
            f"{version}-1",
            "--distribution",
            "unstable",
            "-b",
            "New upstream commits",
        ],
        cwd=tree_root,
    )
    run(["dpkg-source", "-b", str(tree_root)], cwd=source_root)
    return tree_root, dsc_path, version
