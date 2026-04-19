from __future__ import annotations

from . import assets as assets_module
from . import mirrors as mirrors_module
from . import prepare as prepare_module
from . import seed as seed_module
from .constants import (
    SEED_SCHEMA_VERSION,
    base_tgz_path,
    seed_lock_path,
    seed_metadata_path,
)
from .mirrors import pbuilder_mirrorsite, pbuilder_othermirrors
from .seed import base_tgz_is_valid, build_seed_spec, seed_hash

__all__ = [
    "SEED_SCHEMA_VERSION",
    "STAGED_APT_KEYRING_NAME",
    "base_tgz_is_valid",
    "base_tgz_path",
    "build_seed_spec",
    "pbuilder_mirrorsite",
    "pbuilder_othermirrors",
    "prepare_base",
    "prepare_runtime_assets",
    "refresh_feelpp_keyring",
    "resolve_packaging_tree",
    "seed_hash",
    "seed_lock_path",
    "seed_metadata_path",
]

STAGED_APT_KEYRING_NAME = assets_module.STAGED_APT_KEYRING_NAME
prepare_runtime_assets = assets_module.prepare_runtime_assets
refresh_feelpp_keyring = assets_module.refresh_feelpp_keyring
resolve_packaging_tree = assets_module.resolve_packaging_tree


def prepare_base(context: PackagingContext, *, dry_run: bool = False) -> None:
    spec = build_seed_spec(context)
    seed_hash_value = seed_hash(spec)
    base_tgz = base_tgz_path(context)
    metadata_path = seed_metadata_path(context)
    force_refresh = seed_module._force_refresh_requested()
    prepare_runtime_assets(context)

    if dry_run:
        metadata = seed_module._load_seed_metadata(metadata_path)
        base_valid = base_tgz_is_valid(base_tgz)
        if (
            not force_refresh
            and base_valid
            and metadata
            and metadata.get("seed_hash") == seed_hash_value
        ):
            print(f"# reusing seeded pbuilder base {base_tgz} [{seed_hash_value[:12]}]")
            return
        prepare_module._run_prepare_base_shell_with_mirror(
            context,
            mirrorsite=pbuilder_mirrorsite(context),
            othermirrors=pbuilder_othermirrors(context),
            dry_run=True,
        )
        return

    with seed_module._pbuilder_lock(seed_lock_path(context)):
        metadata = seed_module._load_seed_metadata(metadata_path)
        base_valid = base_tgz_is_valid(base_tgz)
        prepare_attempts = seed_module._prepare_max_attempts()
        mirror_candidates = mirrors_module._mirror_site_candidates(context)

        if base_tgz.exists() and not base_valid:
            print(f"--- removing invalid pbuilder base at {base_tgz}")
            seed_module._remove_seed_artifacts(base_tgz, metadata_path)
            metadata = None
            base_valid = False

        if (
            not force_refresh
            and base_valid
            and metadata
            and metadata.get("seed_hash") == seed_hash_value
        ):
            print(
                f"--- reusing seeded pbuilder base at {base_tgz} "
                f"for {context.flavor}/{context.dist} [{seed_hash_value[:12]}]"
            )
            return

        if base_valid and metadata and metadata.get("seed_hash") != seed_hash_value:
            print(
                f"--- invalidating seeded pbuilder base at {base_tgz} "
                f"for {context.flavor}/{context.dist}: seed hash changed"
            )
            seed_module._remove_seed_artifacts(base_tgz, metadata_path)

        last_error: Exception | None = None
        for attempt in range(1, prepare_attempts + 1):
            mirrorsite = mirror_candidates[min(attempt - 1, len(mirror_candidates) - 1)]
            othermirrors = pbuilder_othermirrors(context, mirrorsite=mirrorsite)
            try:
                print(
                    f"--- preparing seeded pbuilder base at {base_tgz} "
                    f"for {context.flavor}/{context.dist} [{seed_hash_value[:12]}] "
                    f"(attempt {attempt}/{prepare_attempts}, mirror {mirrorsite})"
                )
                prepare_module._run_prepare_base_shell_with_mirror(
                    context,
                    mirrorsite=mirrorsite,
                    othermirrors=othermirrors,
                    dry_run=False,
                )
                if not base_tgz_is_valid(base_tgz):
                    raise RuntimeError(
                        f"Prepared pbuilder base archive is missing or invalid: {base_tgz}"
                    )
                seed_module._write_seed_metadata(
                    metadata_path,
                    spec=spec,
                    seed_hash_value=seed_hash_value,
                    base_tgz=base_tgz,
                )
                return
            except subprocess.CalledProcessError as exc:
                last_error = exc
                seed_module._remove_seed_artifacts(base_tgz, metadata_path)
                if attempt == prepare_attempts:
                    raise
                print(
                    f"--- pbuilder base prepare failed for {context.flavor}/{context.dist} "
                    f"(attempt {attempt}/{prepare_attempts}); retrying after cleanup"
                )

        if last_error is not None:
            raise last_error
from ..config import PackagingContext
import subprocess
