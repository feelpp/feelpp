from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess


SOURCE_SUFFIXES = {
    ".c",
    ".cc",
    ".cmake",
    ".cpp",
    ".cxx",
    ".f",
    ".f90",
    ".h",
    ".hh",
    ".hpp",
    ".hxx",
    ".m",
    ".mm",
    ".py",
    ".sh",
}
SOURCE_FILENAMES = {"CMakeLists.txt", "Dockerfile", "Makefile"}
EXCLUDED_PARTS = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".tox",
    ".venv",
    "__pycache__",
    "build",
    "dist",
}
EXCLUDED_PREFIXES = (
    Path("feelpp/contrib"),
    Path("metadata/exports"),
)


@dataclass(frozen=True)
class CodeMetrics:
    approximate_lines_of_code: int
    approximate_bytes: int
    counted_files: int
    method: str = "tracked text source files, excluding generated/build/vendor paths"

    def as_dict(self) -> dict[str, int | str]:
        return {
            "approximate_lines_of_code": self.approximate_lines_of_code,
            "approximate_bytes": self.approximate_bytes,
            "counted_files": self.counted_files,
            "method": self.method,
        }


def _is_relative_to(path: Path, prefix: Path) -> bool:
    try:
        path.relative_to(prefix)
    except ValueError:
        return False
    return True


def _is_source_path(relative_path: Path) -> bool:
    if any(part in EXCLUDED_PARTS for part in relative_path.parts):
        return False
    if any(_is_relative_to(relative_path, prefix) for prefix in EXCLUDED_PREFIXES):
        return False
    return relative_path.name in SOURCE_FILENAMES or relative_path.suffix.lower() in SOURCE_SUFFIXES


def _tracked_paths(repo_root: Path) -> list[Path]:
    try:
        completed = subprocess.run(
            ["git", "ls-files", "-z"],
            cwd=repo_root,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return [
            path.relative_to(repo_root)
            for path in repo_root.rglob("*")
            if path.is_file() and ".git" not in path.relative_to(repo_root).parts
        ]

    names = completed.stdout.decode("utf-8", errors="replace").split("\0")
    return [Path(name) for name in names if name]


def collect_code_metrics(repo_root: str | Path) -> CodeMetrics:
    root = Path(repo_root).expanduser().resolve()
    lines = 0
    bytes_count = 0
    files = 0
    for relative_path in sorted(_tracked_paths(root), key=lambda path: path.as_posix()):
        if not _is_source_path(relative_path):
            continue
        path = root / relative_path
        if not path.is_file():
            continue
        data = path.read_bytes()
        if b"\0" in data:
            continue
        files += 1
        bytes_count += len(data)
        lines += data.count(b"\n")
        if data and not data.endswith(b"\n"):
            lines += 1
    return CodeMetrics(
        approximate_lines_of_code=lines,
        approximate_bytes=bytes_count,
        counted_files=files,
    )
