from __future__ import annotations

import argparse
import json
from pathlib import Path

from ....core.context import discover_repo_root


def _resolve_repo_root(raw: str | None) -> Path:
    return Path(raw).expanduser().resolve() if raw else discover_repo_root()


def _environments_root(repo_root: Path) -> Path:
    return repo_root / "packaging" / "spack" / "environments"


def _environment_manifests(root: Path) -> list[Path]:
    if not root.is_dir():
        raise FileNotFoundError(f"Spack environments directory not found: {root}")
    return sorted(path for path in root.rglob("spack.yaml") if path.is_file())


def _environment_name(root: Path, manifest: Path) -> str:
    return manifest.relative_to(root).parent.as_posix()


def _environment_payload(root: Path, manifest: Path) -> dict[str, str]:
    name = _environment_name(root, manifest)
    return {
        "name": name,
        "environment_root": str(manifest.parent),
        "manifest_path": str(manifest),
        "status": "legacy" if name == "openmpi4" else "supported",
    }


def _find_environment_manifest(root: Path, name: str) -> Path:
    for manifest in _environment_manifests(root):
        if _environment_name(root, manifest) == name:
            return manifest
    raise ValueError(f"Unknown Spack environment: {name}")


def command_spack_env_list(args: argparse.Namespace) -> int:
    repo_root = _resolve_repo_root(args.repo_root)
    root = _environments_root(repo_root)
    manifests = _environment_manifests(root)
    payload = {
        "repo_root": str(repo_root),
        "environments_root": str(root),
        "environments": [_environment_payload(root, manifest) for manifest in manifests],
    }
    print(json.dumps(payload, indent=2))
    return 0


def command_spack_env_show(args: argparse.Namespace) -> int:
    repo_root = _resolve_repo_root(args.repo_root)
    root = _environments_root(repo_root)
    manifest = _find_environment_manifest(root, args.name)
    payload = {
        "repo_root": str(repo_root),
        "environment": _environment_payload(root, manifest),
    }
    print(json.dumps(payload, indent=2))
    return 0


def _add_repo_root_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--repo-root", help="Path to the Feel++ repository root")


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
    env_list_parser.set_defaults(func=command_spack_env_list)

    env_show_parser = env_subparsers.add_parser("show", help="Show one packaged Spack environment")
    _add_repo_root_argument(env_show_parser)
    env_show_parser.add_argument(
        "name",
        help="Environment name relative to packaging/spack/environments",
    )
    env_show_parser.set_defaults(func=command_spack_env_show)
