from __future__ import annotations

from pathlib import Path
import json
import os
import shlex
import shutil
import stat
import subprocess

from ..shell import run_checked
from .constants import (
    CONTAINER_APTLY_ROOT,
    CONTAINER_GNUPG_HOME,
    CONTAINER_STAGED_GNUPG_HOME,
    DEFAULT_APT_SIGNING_KEY,
    DEFAULT_STAGED_APT_KEYRING,
)


def _host_aptly_config_path() -> Path | None:
    candidates: list[Path] = []
    explicit = os.getenv("FEELPP_APTLY_CONFIG")
    if explicit:
        candidates.append(Path(explicit).expanduser())
    candidates.extend(
        [
            Path.home() / ".aptly.conf",
            Path.home() / ".config" / "aptly.conf",
            Path.home() / ".config" / "aptly" / "aptly.conf",
        ]
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return None


def _stage_aptly_config(host_config_path: Path, *, job_root: Path) -> Path:
    staged_config = job_root / "aptly.conf"
    data = json.loads(host_config_path.read_text(encoding="utf-8"))
    data["rootDir"] = str(CONTAINER_APTLY_ROOT)
    staged_config.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    staged_config.chmod(0o600)
    return staged_config


def _host_gnupg_home() -> Path | None:
    candidate = Path(os.getenv("GNUPGHOME", str(Path.home() / ".gnupg"))).expanduser()
    if not candidate.is_dir():
        return None
    return candidate.resolve()


def _stage_host_gnupg_home(host_gnupg_home: Path, *, job_root: Path) -> Path:
    staged_home = job_root / "gnupg-home"
    if staged_home.exists():
        shutil.rmtree(staged_home)
    staged_home.mkdir(parents=True, exist_ok=True)
    staged_home.chmod(0o700)

    for current_root, dirnames, filenames in os.walk(host_gnupg_home):
        root_path = Path(current_root)
        relative_root = root_path.relative_to(host_gnupg_home)
        target_root = staged_home / relative_root
        target_root.mkdir(parents=True, exist_ok=True)
        target_root.chmod(0o700)

        filtered_dirnames: list[str] = []
        for dirname in dirnames:
            source_dir = root_path / dirname
            try:
                mode = source_dir.lstat().st_mode
            except FileNotFoundError:
                continue
            if stat.S_ISSOCK(mode) or source_dir.is_symlink() or dirname.startswith(".#lk"):
                continue
            filtered_dirnames.append(dirname)
        dirnames[:] = filtered_dirnames

        for filename in filenames:
            source = root_path / filename
            try:
                mode = source.lstat().st_mode
            except FileNotFoundError:
                continue
            if stat.S_ISSOCK(mode) or source.is_symlink() or filename.startswith("S.gpg-agent") or filename.startswith(".#lk"):
                continue
            destination = target_root / filename
            shutil.copy2(source, destination)
            destination.chmod(0o600)

    return staged_home


def _feelpp_apt_signing_key() -> str:
    return os.getenv("FEELPP_APT_GPG_KEY") or os.getenv("GPG_KEY") or DEFAULT_APT_SIGNING_KEY


def _stage_host_feelpp_apt_key(*, job_root: Path) -> Path | None:
    staged_key = job_root / DEFAULT_STAGED_APT_KEYRING
    explicit_key_file = os.getenv("FEELPP_APT_KEY_FILE")
    if explicit_key_file:
        source = Path(explicit_key_file).expanduser().resolve()
        if not source.is_file():
            return None
        shutil.copyfile(source, staged_key)
        staged_key.chmod(0o644)
        return staged_key

    try:
        run_checked(
            [
                "gpg",
                "--batch",
                "--yes",
                "--output",
                str(staged_key),
                "--export",
                _feelpp_apt_signing_key(),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None

    if not staged_key.is_file():
        return None
    staged_key.chmod(0o644)
    return staged_key


def _gpg_bootstrap_commands(*, host_gnupg_mounted: bool) -> list[str]:
    commands = [
        f"rm -rf {shlex.quote(str(CONTAINER_GNUPG_HOME))}",
        f"install -d -m 700 {shlex.quote(str(CONTAINER_GNUPG_HOME))}",
        f"export GNUPGHOME={shlex.quote(str(CONTAINER_GNUPG_HOME))}",
    ]
    if not host_gnupg_mounted:
        return commands

    commands.extend(
        [
            f"if [ -d {shlex.quote(str(CONTAINER_STAGED_GNUPG_HOME))} ]; then",
            f"  tar -C {shlex.quote(str(CONTAINER_STAGED_GNUPG_HOME))} -cf - . | tar -C {shlex.quote(str(CONTAINER_GNUPG_HOME))} -xf -",
            f"  chown -R root:root {shlex.quote(str(CONTAINER_GNUPG_HOME))}",
            f"  if find {shlex.quote(str(CONTAINER_STAGED_GNUPG_HOME))} -maxdepth 1 -type f \\( -name '*.asc' -o -name '*.gpg' \\) -print -quit | grep -q .; then",
            f"    find {shlex.quote(str(CONTAINER_STAGED_GNUPG_HOME))} -maxdepth 1 -type f \\( -name '*.asc' -o -name '*.gpg' \\) -print0 | xargs -0 -r gpg --batch --import >/dev/null 2>&1 || true",
            "  fi",
            "fi",
            f"find {shlex.quote(str(CONTAINER_GNUPG_HOME))} -type d -exec chmod 700 {{}} +",
            f"find {shlex.quote(str(CONTAINER_GNUPG_HOME))} -type f -exec chmod u+rw,go-rwx {{}} +",
        ]
    )
    return commands
