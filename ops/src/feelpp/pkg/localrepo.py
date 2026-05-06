from __future__ import annotations

from pathlib import Path
import gzip
import shutil
import subprocess

from .config import PackagingContext
from .workspace import ensure_workspace


def stage_outputs(context: PackagingContext, result_dir: Path) -> dict[str, int]:
    ensure_workspace(context)
    copied_files = 0
    repo_packages = 0

    for path in sorted(result_dir.iterdir()):
        if not path.is_file():
            continue
        shutil.copy2(path, context.artifacts_dir / path.name)
        copied_files += 1
        if path.suffix in {".deb", ".udeb"}:
            package_name = path.name.split("_", 1)[0]
            for existing in context.local_repo_dir.glob(f"{package_name}_*{path.suffix}"):
                if existing.name != path.name:
                    existing.unlink()
            shutil.copy2(path, context.local_repo_dir / path.name)
            repo_packages += 1

    if repo_packages:
        packages = subprocess.check_output(
            ["dpkg-scanpackages", ".", "/dev/null"],
            cwd=context.local_repo_dir,
            text=True,
        )
        packages_path = context.local_repo_dir / "Packages"
        packages_path.write_text(packages, encoding="utf-8")
        with gzip.open(context.local_repo_dir / "Packages.gz", "wt", encoding="utf-8") as handle:
            handle.write(packages)

    return {"copied_files": copied_files, "repo_packages": repo_packages}
