from __future__ import annotations

from pathlib import Path

from ..config import PackagingContext


SEED_SCHEMA_VERSION = 1
STAGED_APT_KEYRING_NAME = "feelpp-archive-keyring.gpg"


def base_tgz_path(context: PackagingContext) -> Path:
    return context.pbuilder_root / f"{context.dist}-base.tgz"


def seed_metadata_path(context: PackagingContext) -> Path:
    return context.pbuilder_root / f"{context.dist}-base.seed.json"


def seed_lock_path(context: PackagingContext) -> Path:
    return context.pbuilder_root / f"{context.dist}-base.lock"
