from __future__ import annotations

import contextlib
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Iterator

from ..config import PackagingContext
from ..debian import collect_seed_build_dependencies, control_paths_for_context
from ..shell import run_probe
from .constants import SEED_SCHEMA_VERSION
from .mirrors import _mirror_site, _mirror_site_candidates


TRUTHY = {"1", "true", "yes", "on"}
DEFAULT_PREPARE_MAX_ATTEMPTS = 4
DEFAULT_BUILDDEPS_BATCH_SIZE = 8


def _force_refresh_requested() -> bool:
    return os.getenv("FEELPP_PKG_PBUILDER_REFRESH", "").strip().lower() in TRUTHY


def _seed_input_paths(context: PackagingContext) -> list[Path]:
    package_dir = Path(__file__).resolve().parent
    generator_inputs = [
        package_dir / "__init__.py",
        package_dir / "assets.py",
        package_dir / "constants.py",
        package_dir / "mirrors.py",
        package_dir / "prepare.py",
        package_dir / "seed.py",
        Path(__file__).resolve().parents[1] / "debian.py",
    ]
    paths = [context.pbuilder_config]
    top_level_hooks = sorted(
        hook
        for hook in context.pbuilder_source_hookdir.glob("[A-Z0-9]*-*")
        if hook.name != "E10-feelpp-builddeps"
    )
    keyrings = sorted((context.pbuilder_source_hookdir / "keyrings").glob("*.gpg.b64"))
    controls = control_paths_for_context(context)
    manifest = [context.manifest_path]
    return [*paths, *top_level_hooks, *keyrings, *manifest, *controls, *generator_inputs]


def build_seed_spec(context: PackagingContext) -> dict[str, object]:
    inputs: dict[str, str] = {}
    for path in _seed_input_paths(context):
        if not path.is_file():
            continue
        try:
            key = path.relative_to(context.repo_root).as_posix()
        except ValueError:
            key = f"external::{path}"
        inputs[key] = hashlib.sha256(path.read_bytes()).hexdigest()

    return {
        "schema_version": SEED_SCHEMA_VERSION,
        "dist": context.dist,
        "flavor": context.flavor,
        "channel": context.channel,
        "mirror_site": _mirror_site(context),
        "mirror_candidates": _mirror_site_candidates(context),
        "inputs": inputs,
    }


def seed_hash(spec: dict[str, object]) -> str:
    payload = json.dumps(spec, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _load_builddeps(context: PackagingContext) -> list[str]:
    entries = collect_seed_build_dependencies(context)
    if not entries:
        raise ValueError(f"Derived build dependency list is empty for {context.dist}")
    return entries


def _builddeps_batch_size() -> int:
    raw = os.getenv("FEELPP_PBUILDER_BUILDDEPS_BATCH_SIZE", "").strip()
    if not raw:
        return DEFAULT_BUILDDEPS_BATCH_SIZE
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"Invalid FEELPP_PBUILDER_BUILDDEPS_BATCH_SIZE value: {raw}") from exc
    if value < 1:
        raise ValueError("FEELPP_PBUILDER_BUILDDEPS_BATCH_SIZE must be >= 1")
    return value


def _prepare_max_attempts() -> int:
    raw = os.getenv("FEELPP_PBUILDER_PREPARE_MAX_ATTEMPTS", "").strip()
    if not raw:
        return DEFAULT_PREPARE_MAX_ATTEMPTS
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"Invalid FEELPP_PBUILDER_PREPARE_MAX_ATTEMPTS value: {raw}") from exc
    if value < 1:
        raise ValueError("FEELPP_PBUILDER_PREPARE_MAX_ATTEMPTS must be >= 1")
    return value


def _batched(items: list[str], size: int) -> list[list[str]]:
    return [items[index : index + size] for index in range(0, len(items), size)]


def _load_seed_metadata(path: Path) -> dict[str, object] | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _write_seed_metadata(
    path: Path,
    *,
    spec: dict[str, object],
    seed_hash_value: str,
    base_tgz: Path,
) -> None:
    payload = {
        "schema_version": SEED_SCHEMA_VERSION,
        "seed_hash": seed_hash_value,
        "base_tgz": str(base_tgz),
        "prepared_at": datetime.now(timezone.utc).isoformat(),
        "spec": spec,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def base_tgz_is_valid(path: Path) -> bool:
    if not path.is_file():
        return False
    return run_probe(
        ["tar", "-tzf", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def _remove_seed_artifacts(base_tgz: Path, metadata_path: Path) -> None:
    if base_tgz.exists():
        base_tgz.unlink()
    if metadata_path.exists():
        metadata_path.unlink()


@contextlib.contextmanager
def _pbuilder_lock(path: Path) -> Iterator[None]:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
