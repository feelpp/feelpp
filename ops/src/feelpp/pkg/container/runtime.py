from __future__ import annotations

from pathlib import Path
import os
import tempfile

from ..config import PackagingContext
from .constants import DEFAULT_DOCKER_IMAGE, DEFAULT_DOCKER_IMAGES_BY_CONTEXT


def default_image_for_context(context: PackagingContext) -> str:
    override = os.getenv("FEELPP_PKG_DOCKER_IMAGE")
    if override:
        return override
    return DEFAULT_DOCKER_IMAGES_BY_CONTEXT.get((context.flavor, context.dist), DEFAULT_DOCKER_IMAGE)


def _preferred_site_cache_home() -> Path | None:
    site_root = Path("/nvme0/cemosis")
    if site_root.is_dir() and os.access(site_root, os.W_OK):
        return site_root / ".cache"
    return None


def _default_cache_home() -> Path:
    preferred_site_cache_home = _preferred_site_cache_home()
    if preferred_site_cache_home is not None:
        return preferred_site_cache_home

    xdg_cache_home = os.getenv("XDG_CACHE_HOME")
    if xdg_cache_home:
        return Path(xdg_cache_home).expanduser()

    home = os.getenv("HOME")
    if home:
        return Path(home).expanduser() / ".cache"

    return Path(tempfile.gettempdir()) / "feelpp-cache"


def default_state_root() -> Path:
    explicit = os.getenv("FEELPP_PKG_DOCKER_STATE_ROOT")
    if explicit:
        return Path(explicit).expanduser().resolve()

    explicit_base = os.getenv("FEELPP_PKG_DOCKER_STATE_ROOT_BASE")
    if explicit_base:
        return Path(explicit_base).expanduser().resolve()

    return (_default_cache_home() / "feelpp-pkg" / "docker").resolve()


def _bootstrap_apt_packages(context: PackagingContext) -> list[str]:
    packages = ["pkgconf", "arch-test"]
    if context.flavor == "debian":
        packages.append("debian-archive-keyring")
    return packages
