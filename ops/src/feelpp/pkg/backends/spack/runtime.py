from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shlex
import shutil
import tempfile

from ...core.context import discover_repo_root


def resolve_repo_root(raw: str | None = None) -> Path:
    return Path(raw).expanduser().resolve() if raw else discover_repo_root()


def spack_metadata_root(repo_root: Path, raw: str | None = None) -> Path:
    explicit = raw or os.getenv("FEELPP_SPACK_METADATA_ROOT")
    if explicit:
        return Path(explicit).expanduser().resolve()

    candidates = (
        repo_root / "packaging" / "spack",
        repo_root / "ops" / "packaging" / "spack",
    )
    for candidate in candidates:
        if (candidate / "environments").is_dir():
            return candidate
    return candidates[0]


def environments_root(repo_root: Path, raw_metadata_root: str | None = None) -> Path:
    return spack_metadata_root(repo_root, raw_metadata_root) / "environments"


def _repo_key(repo_root: Path) -> str:
    digest = hashlib.sha1(str(repo_root).encode("utf-8")).hexdigest()[:8]
    return f"{repo_root.name}-{digest}"


def _default_cache_home() -> Path:
    xdg_cache_home = os.getenv("XDG_CACHE_HOME")
    if xdg_cache_home:
        return Path(xdg_cache_home).expanduser()

    home = os.getenv("HOME")
    if home:
        return Path(home).expanduser() / ".cache"

    return Path(tempfile.gettempdir()) / "feelpp-cache"


def _default_config_home() -> Path:
    xdg_config_home = os.getenv("XDG_CONFIG_HOME")
    if xdg_config_home:
        return Path(xdg_config_home).expanduser()

    home = os.getenv("HOME")
    if home:
        return Path(home).expanduser() / ".config"

    return Path(tempfile.gettempdir()) / "feelpp-config"


def default_user_cache_path(repo_root: Path) -> Path:
    return (_default_cache_home() / "feelpp-spack" / _repo_key(repo_root)).resolve()


def default_user_config_path(repo_root: Path) -> Path:
    return (_default_config_home() / "feelpp-spack" / _repo_key(repo_root)).resolve()


def resolve_user_cache_path(repo_root: Path, raw: str | None = None) -> Path:
    explicit = raw or os.getenv("SPACK_USER_CACHE_PATH") or os.getenv("FEELPP_SPACK_USER_CACHE_PATH")
    if explicit:
        return Path(explicit).expanduser().resolve()
    return default_user_cache_path(repo_root)


def resolve_user_config_path(repo_root: Path, raw: str | None = None) -> Path:
    explicit = raw or os.getenv("SPACK_USER_CONFIG_PATH") or os.getenv("FEELPP_SPACK_USER_CONFIG_PATH")
    if explicit:
        return Path(explicit).expanduser().resolve()
    return default_user_config_path(repo_root)


def resolve_spack_root(repo_root: Path, raw: str | None = None) -> Path:
    explicit = raw or os.getenv("FEELPP_SPACK_ROOT") or os.getenv("SPACK_ROOT")
    if explicit:
        return Path(explicit).expanduser().resolve()
    return (resolve_user_cache_path(repo_root) / "spack").resolve()


def spack_executable(spack_root: Path | None = None) -> Path | None:
    if spack_root is not None:
        candidate = spack_root / "bin" / "spack"
        if candidate.is_file():
            return candidate

    found = shutil.which("spack")
    if found:
        return Path(found).resolve()

    if spack_root is not None:
        return spack_root / "bin" / "spack"
    return None


def shell_quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def spack_env_overrides(
    *,
    user_config_path: Path,
    user_cache_path: Path,
) -> dict[str, str]:
    return {
        "SPACK_USER_CONFIG_PATH": str(user_config_path),
        "SPACK_USER_CACHE_PATH": str(user_cache_path),
    }
