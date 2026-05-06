from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import yaml

from ..runtime import environments_root, resolve_repo_root, resolve_user_cache_path, spack_metadata_root


def _resolve_repo_root(raw: str | None) -> Path:
    return resolve_repo_root(raw)


def _environments_root(repo_root: Path) -> Path:
    return environments_root(repo_root)


def _environment_manifests(root: Path) -> list[Path]:
    if not root.is_dir():
        raise FileNotFoundError(f"Spack environments directory not found: {root}")
    return sorted(path for path in root.rglob("spack.yaml") if path.is_file())


def _environment_name(root: Path, manifest: Path) -> str:
    return manifest.relative_to(root).parent.as_posix()


def _environment_status(name: str) -> str:
    return "legacy" if name == "openmpi4" else "supported"


def _environment_aliases(root: Path, manifest: Path) -> list[str]:
    name = _environment_name(root, manifest)
    aliases: list[str] = []
    leaf_name = Path(name).name
    if leaf_name != name:
        aliases.append(leaf_name)
    if name.startswith("cpu/"):
        cpu_name = name.removeprefix("cpu/")
        if cpu_name not in aliases:
            aliases.append(cpu_name)
    return aliases


def _load_manifest_payload(manifest: Path) -> dict:
    payload = yaml.safe_load(manifest.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"Spack manifest is not a YAML mapping: {manifest}")
    return payload


def _expand_view_path(value: str, *, manifest: Path, user_cache_path: Path) -> Path:
    expanded = (
        value.replace("$user_cache_path", str(user_cache_path))
        .replace("${user_cache_path}", str(user_cache_path))
        .replace("$env", str(manifest.parent))
        .replace("${env}", str(manifest.parent))
    )
    return Path(os.path.expandvars(expanded)).expanduser().resolve()


def _view_entry(
    *,
    name: str,
    spec: object,
    manifest: Path,
    user_cache_path: Path,
) -> dict[str, object]:
    path: Path | None = None
    if spec is True:
        path = manifest.parent / ".spack-env" / "view"
    elif isinstance(spec, str):
        path = _expand_view_path(spec, manifest=manifest, user_cache_path=user_cache_path)
    elif isinstance(spec, dict):
        root = spec.get("root") or spec.get("path")
        if isinstance(root, str):
            path = _expand_view_path(root, manifest=manifest, user_cache_path=user_cache_path)

    return {
        "name": name,
        "spec": spec,
        "path": str(path) if path is not None else None,
        "exists": path.exists() if path is not None else None,
    }


def _environment_views(manifest: Path, payload: dict, user_cache_path: Path) -> list[dict[str, object]]:
    spack_payload = payload.get("spack", {})
    if not isinstance(spack_payload, dict) or "view" not in spack_payload:
        return []

    view_spec = spack_payload["view"]
    if view_spec in (False, None):
        return []
    if isinstance(view_spec, dict):
        return [
            _view_entry(
                name=str(name),
                spec=spec,
                manifest=manifest,
                user_cache_path=user_cache_path,
            )
            for name, spec in view_spec.items()
        ]
    return [
        _view_entry(
            name="default",
            spec=view_spec,
            manifest=manifest,
            user_cache_path=user_cache_path,
        )
    ]


def _repo_root_from_environments_root(root: Path) -> Path:
    for candidate in root.parents:
        if root == candidate / "packaging" / "spack" / "environments":
            return candidate
        if root == candidate / "ops" / "packaging" / "spack" / "environments":
            return candidate
    return root.parents[2]


def _environment_payload(
    root: Path,
    manifest: Path,
    *,
    user_cache_path: Path | None = None,
) -> dict[str, object]:
    name = _environment_name(root, manifest)
    resolved_user_cache_path = user_cache_path or resolve_user_cache_path(
        _repo_root_from_environments_root(root)
    )
    payload: dict[str, object] = {
        "name": name,
        "aliases": _environment_aliases(root, manifest),
        "environment_root": str(manifest.parent),
        "manifest_path": str(manifest),
        "status": _environment_status(name),
    }
    try:
        manifest_payload = _load_manifest_payload(manifest)
    except Exception as exc:  # noqa: BLE001
        payload["manifest_error"] = str(exc)
        payload["views"] = []
        return payload
    payload["views"] = _environment_views(manifest, manifest_payload, resolved_user_cache_path)
    return {
        **payload,
    }


def _find_environment_manifest(root: Path, name: str) -> Path:
    normalized_name = name.strip("/")
    by_name = {
        _environment_name(root, manifest): manifest
        for manifest in _environment_manifests(root)
    }
    if normalized_name in by_name:
        return by_name[normalized_name]

    candidates = [
        (environment_name, manifest)
        for environment_name, manifest in by_name.items()
        if Path(environment_name).name == normalized_name
        or (
            environment_name.startswith("cpu/")
            and environment_name.removeprefix("cpu/") == normalized_name
        )
    ]
    if len(candidates) == 1:
        return candidates[0][1]
    if len(candidates) > 1:
        names = ", ".join(environment_name for environment_name, _ in candidates)
        raise ValueError(f"Ambiguous Spack environment {name!r}; matches: {names}")

    for manifest in _environment_manifests(root):
        if _environment_name(root, manifest) == normalized_name:
            return manifest
    raise ValueError(f"Unknown Spack environment: {name}")


def command_spack_env_list(args: argparse.Namespace) -> int:
    repo_root = _resolve_repo_root(args.repo_root)
    root = _environments_root(repo_root)
    metadata_root = spack_metadata_root(repo_root)
    user_cache_path = resolve_user_cache_path(repo_root, getattr(args, "user_cache_path", None))
    manifests = _environment_manifests(root)
    payload = {
        "repo_root": str(repo_root),
        "spack_metadata_root": str(metadata_root),
        "environments_root": str(root),
        "spack_user_cache_path": str(user_cache_path),
        "environments": [
            _environment_payload(root, manifest, user_cache_path=user_cache_path)
            for manifest in manifests
        ],
    }
    print(json.dumps(payload, indent=2))
    return 0


def command_spack_env_show(args: argparse.Namespace) -> int:
    repo_root = _resolve_repo_root(args.repo_root)
    root = _environments_root(repo_root)
    manifest = _find_environment_manifest(root, args.name)
    user_cache_path = resolve_user_cache_path(repo_root, getattr(args, "user_cache_path", None))
    payload = {
        "repo_root": str(repo_root),
        "spack_metadata_root": str(spack_metadata_root(repo_root)),
        "environment": _environment_payload(root, manifest, user_cache_path=user_cache_path),
    }
    print(json.dumps(payload, indent=2))
    return 0


def _add_repo_root_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--repo-root", help="Path to the Feel++ repository root")


def _add_user_cache_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--user-cache-path",
        help="SPACK_USER_CACHE_PATH used to resolve environment view paths",
    )


def register_env_commands(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    env_parser = subparsers.add_parser(
        "env",
        help="Inspect repository-owned Spack environments",
    )
    env_subparsers = env_parser.add_subparsers(dest="spack_env_command", required=True)

    env_list_parser = env_subparsers.add_parser("list", help="List packaged Spack environments")
    _add_repo_root_argument(env_list_parser)
    _add_user_cache_argument(env_list_parser)
    env_list_parser.set_defaults(func=command_spack_env_list)

    env_show_parser = env_subparsers.add_parser("show", help="Show one packaged Spack environment")
    _add_repo_root_argument(env_show_parser)
    _add_user_cache_argument(env_show_parser)
    env_show_parser.add_argument(
        "name",
        help="Environment name relative to packaging/spack/environments",
    )
    env_show_parser.set_defaults(func=command_spack_env_show)
