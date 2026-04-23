from __future__ import annotations

from pathlib import Path
from typing import Any

from feelpp.pkg.config import discover_repo_root
from feelpp.ops.version.repository import VersionRepository

from . import PUBLIC_METADATA_RELATIVE_PATH
from .metrics import CodeMetrics
from .schema import validate_public_metadata
from .yamlio import load_yaml, write_yaml


class PublicMetadataRepository:
    def __init__(
        self,
        repo_root: str | Path | None = None,
        metadata_path: str | Path | None = None,
    ) -> None:
        self.repo_root = Path(repo_root).expanduser().resolve() if repo_root else discover_repo_root()
        self.metadata_path = (
            Path(metadata_path).expanduser().resolve()
            if metadata_path
            else self.repo_root / PUBLIC_METADATA_RELATIVE_PATH
        )

    def relative_path(self, path: Path) -> str:
        try:
            return path.resolve().relative_to(self.repo_root).as_posix()
        except ValueError:
            return path.as_posix()

    def read_public_metadata(self) -> dict[str, Any]:
        payload = load_yaml(self.metadata_path)
        validate_public_metadata(payload)
        return payload

    def write_public_metadata(self, payload: dict[str, Any]) -> None:
        validate_public_metadata(payload)
        write_yaml(self.metadata_path, payload)

    def write_metrics(self, metrics: CodeMetrics) -> dict[str, Any]:
        payload = self.read_public_metadata()
        updated_metrics = dict(payload.get("metrics") or {})
        updated_metrics.update(metrics.as_dict())
        payload["metrics"] = updated_metrics
        self.write_public_metadata(payload)
        return payload

    def release_metadata(self) -> dict[str, Any]:
        state = VersionRepository(repo_root=self.repo_root).read_state()
        return {
            "version": str(state.canonical_upstream_version()),
            "cmake_versions": [
                {
                    "name": record.name,
                    "path": self.relative_path(record.path),
                    "version": str(record.version),
                }
                for record in sorted(state.cmake_versions, key=lambda item: item.name)
            ],
            "public_metadata_versions": [
                {
                    "name": record.name,
                    "path": self.relative_path(record.path),
                    "version": str(record.version),
                }
                for record in sorted(state.metadata_versions, key=lambda item: item.name)
            ],
            "package_versions": [
                {
                    "component": record.component,
                    "dist": record.dist,
                    "version": str(record.version),
                    "source_name": record.source_name,
                }
                for record in sorted(
                    state.package_versions,
                    key=lambda item: (item.component, item.dist, str(item.version)),
                )
            ],
        }
