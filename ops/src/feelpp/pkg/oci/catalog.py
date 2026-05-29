from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path

from ..core.context import discover_repo_root


IMAGE_PROFILE = "images"


def _as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


@dataclass(frozen=True)
class ImageTarget:
    target: str
    flavor: str
    dist: str
    version: str
    docker: bool
    continue_on_error: bool
    image_backend: str
    image_strategy: str
    base_image: str
    oci_dist: str
    component_set: str | None = None
    spack_environment: str | None = None
    metadata: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_row(cls, target: str, row: dict[str, object]) -> "ImageTarget":
        normalized = {str(key): str(value) for key, value in row.items()}
        known_keys = {
            "flavor",
            "dist",
            "version",
            "docker",
            "continue-on-error",
            "image_backend",
            "image_strategy",
            "base_image",
            "oci_dist",
            "component_set",
            "spack_environment",
        }
        extra = {
            key: value
            for key, value in normalized.items()
            if key not in known_keys
        }
        return cls(
            target=str(target),
            flavor=normalized.get("flavor", ""),
            dist=normalized.get("dist", ""),
            version=normalized.get("version", ""),
            docker=_as_bool(normalized.get("docker", "false")),
            continue_on_error=_as_bool(normalized.get("continue-on-error", "false")),
            image_backend=normalized.get("image_backend", ""),
            image_strategy=normalized.get("image_strategy", ""),
            base_image=normalized.get("base_image", ""),
            oci_dist=normalized.get("oci_dist", ""),
            component_set=normalized.get("component_set") or None,
            spack_environment=normalized.get("spack_environment") or None,
            metadata=extra,
        )

    def as_dict(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "target": self.target,
            "flavor": self.flavor,
            "dist": self.dist,
            "version": self.version,
            "docker": self.docker,
            "continue-on-error": self.continue_on_error,
            "image_backend": self.image_backend,
            "image_strategy": self.image_strategy,
            "base_image": self.base_image,
            "oci_dist": self.oci_dist,
        }
        if self.component_set:
            payload["component_set"] = self.component_set
        if self.spack_environment:
            payload["spack_environment"] = self.spack_environment
        payload.update(self.metadata)
        return payload


def _resolve_repo_root(raw: str | None) -> Path:
    return Path(raw).expanduser().resolve() if raw else discover_repo_root()


def plan_ci_path(repo_root: Path) -> Path:
    return repo_root / ".github" / "plan-ci.json"


def load_image_profile_config(repo_root: Path, *, profile: str = IMAGE_PROFILE) -> dict[str, object]:
    payload = json.loads(plan_ci_path(repo_root).read_text(encoding="utf-8"))
    profiles = payload.get("profiles", {})
    if not isinstance(profiles, dict):
        raise ValueError(f"Invalid profiles payload in {plan_ci_path(repo_root)}")
    profile_payload = profiles.get(profile, {})
    if not isinstance(profile_payload, dict):
        raise ValueError(f"Invalid profile '{profile}' in {plan_ci_path(repo_root)}")
    return profile_payload


def load_image_catalog(
    repo_root: Path,
    *,
    profile: str = IMAGE_PROFILE,
) -> dict[str, ImageTarget]:
    profile_payload = load_image_profile_config(repo_root, profile=profile)
    catalog = profile_payload.get("catalog", {})
    if not isinstance(catalog, dict):
        raise ValueError(f"Invalid catalog in {plan_ci_path(repo_root)} for profile '{profile}'")
    return {
        str(target): ImageTarget.from_row(str(target), row)
        for target, row in catalog.items()
        if isinstance(row, dict)
    }


def list_image_targets(
    repo_root: Path | str | None = None,
    *,
    profile: str = IMAGE_PROFILE,
    backend: str | None = None,
) -> list[ImageTarget]:
    resolved_repo_root = _resolve_repo_root(str(repo_root) if repo_root else None)
    catalog = load_image_catalog(resolved_repo_root, profile=profile)
    rows = sorted(catalog.values(), key=lambda item: item.target)
    if backend:
        rows = [row for row in rows if row.image_backend == backend]
    return rows


def get_image_target(
    target: str,
    *,
    repo_root: Path | str | None = None,
    profile: str = IMAGE_PROFILE,
    backend: str | None = None,
) -> ImageTarget:
    resolved_repo_root = _resolve_repo_root(str(repo_root) if repo_root else None)
    catalog = load_image_catalog(resolved_repo_root, profile=profile)
    try:
        row = catalog[target]
    except KeyError as exc:
        known = ", ".join(sorted(catalog))
        raise ValueError(f"Unknown image target: {target}. Known targets: {known}") from exc
    if backend and row.image_backend != backend:
        raise ValueError(
            f"Image target {target} uses backend '{row.image_backend}', expected '{backend}'"
        )
    return row
