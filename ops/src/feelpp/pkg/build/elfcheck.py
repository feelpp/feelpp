from __future__ import annotations

from pathlib import Path
import re
import subprocess
import tempfile


_NEEDED_RE = re.compile(r"Shared library: \[(?P<name>[^\]]+)\]")


def _is_dev_package(path: Path) -> bool:
    package_name = path.name.split("_", 1)[0]
    return package_name.endswith("-dev")


def _is_elf(path: Path) -> bool:
    if not path.is_file() or path.is_symlink():
        return False
    try:
        return path.read_bytes().startswith(b"\x7fELF")
    except OSError:
        return False


def _invalid_needed_entries(output: str) -> list[str]:
    invalid: list[str] = []
    for match in _NEEDED_RE.finditer(output):
        name = match.group("name")
        if name.startswith("libfeelpp_") and name.endswith(".so"):
            invalid.append(name)
    return invalid


def validate_runtime_linkage(result_dir: Path) -> None:
    for package in sorted(result_dir.glob("*.deb")):
        if _is_dev_package(package):
            continue

        with tempfile.TemporaryDirectory(prefix="feelpp-pkg-elf-") as tmpdir:
            extract_root = Path(tmpdir)
            subprocess.run(
                ["dpkg-deb", "-x", str(package), str(extract_root)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            for candidate in sorted(extract_root.rglob("*")):
                if not _is_elf(candidate):
                    continue
                completed = subprocess.run(
                    ["readelf", "-d", str(candidate)],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                invalid = _invalid_needed_entries(completed.stdout)
                if invalid:
                    raise RuntimeError(
                        "Runtime package contains unversioned Feel++ linkage: "
                        + ", ".join(invalid)
                    )
