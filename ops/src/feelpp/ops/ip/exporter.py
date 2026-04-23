from __future__ import annotations

from pathlib import Path
from typing import Any

from .repository import PublicMetadataRepository
from .yamlio import dump_yaml


EXPORT_SCHEMA_VERSION = 1
EXPORT_KIND = "feelpp-public-app-bundle"


def build_public_bundle(repository: PublicMetadataRepository) -> dict[str, Any]:
    metadata = repository.read_public_metadata()
    return {
        "schema_version": EXPORT_SCHEMA_VERSION,
        "kind": EXPORT_KIND,
        "source": {
            "metadata_path": repository.relative_path(repository.metadata_path),
            "metadata_schema_version": metadata["schema_version"],
        },
        "software": metadata["software"],
        "repository": metadata["repository"],
        "languages": metadata["languages"],
        "build_tools": metadata["build_tools"],
        "public_license": metadata["public_license"],
        "architecture": metadata["architecture"],
        "public_funding": metadata["public_funding"],
        "public_outputs": metadata["public_outputs"],
        "metrics": metadata["metrics"],
        "release": repository.release_metadata(),
        "private_boundary": metadata["private_boundary"],
    }


def render_public_bundle(payload: dict[str, Any], export_format: str) -> str:
    if export_format == "json":
        import json

        return json.dumps(payload, indent=2, sort_keys=False) + "\n"
    if export_format == "yaml":
        return dump_yaml(payload)
    raise ValueError(f"Unsupported export format: {export_format}")


def write_public_bundle(path: Path, payload: dict[str, Any], export_format: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_public_bundle(payload, export_format), encoding="utf-8")
