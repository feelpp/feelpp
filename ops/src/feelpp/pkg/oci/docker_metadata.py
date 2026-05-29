from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class DockerDistributionVersion:
    family: str
    name: str
    version: str
    base_image: str
    template: str
    platforms: tuple[str, ...]
    experimental: bool
    openmpi_packages: tuple[str, ...]


@dataclass(frozen=True)
class DockerVariant:
    name: str
    description: str
    layers: tuple[str, ...]
    enable_openmpi: bool
    extra_packages: tuple[str, ...]


@dataclass(frozen=True)
class DockerMetadata:
    root: Path
    distributions: dict[str, DockerDistributionVersion]
    package_groups: dict[str, tuple[str, ...]]
    layers: dict[str, tuple[str, ...]]
    variants: dict[str, DockerVariant]


def _read_yaml(path: Path) -> dict[str, object]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid YAML document: {path}")
    return payload


def docker_metadata_root(repo_root: Path) -> Path:
    return repo_root / "packaging" / "docker"


def load_docker_metadata(repo_root: Path) -> DockerMetadata:
    root = docker_metadata_root(repo_root)
    distributions_payload = _read_yaml(root / "config" / "distributions.yaml")
    packages_payload = _read_yaml(root / "config" / "packages.yaml")
    variants_payload = _read_yaml(root / "config" / "variants.yaml")

    distributions: dict[str, DockerDistributionVersion] = {}
    for family_name, family_payload in (distributions_payload.get("distributions", {}) or {}).items():
        if not isinstance(family_payload, dict):
            continue
        template = str(family_payload.get("template", "")).strip()
        for version_payload in family_payload.get("versions", []) or []:
            if not isinstance(version_payload, dict):
                continue
            version_name = str(version_payload.get("name", "")).strip()
            if not version_name:
                continue
            distributions[f"{family_name}:{version_name}"] = DockerDistributionVersion(
                family=str(family_name),
                name=version_name,
                version=str(version_payload.get("version", "")),
                base_image=str(version_payload.get("base_image", "")),
                template=template,
                platforms=tuple(str(item) for item in version_payload.get("platforms", []) or []),
                experimental=bool(version_payload.get("experimental", False)),
                openmpi_packages=tuple(
                    str(item) for item in version_payload.get("openmpi_packages", []) or []
                ),
            )

    package_groups = {
        str(name): tuple(str(item) for item in items or [])
        for name, items in (packages_payload.get("groups", {}) or {}).items()
        if isinstance(items, list)
    }
    layers = {
        str(name): tuple(str(item) for item in items or [])
        for name, items in (packages_payload.get("layers", {}) or {}).items()
        if isinstance(items, list)
    }
    variants = {
        str(name): DockerVariant(
            name=str(name),
            description=str(payload.get("description", "")),
            layers=tuple(str(item) for item in payload.get("layers", []) or []),
            enable_openmpi=bool(payload.get("enable_openmpi", True)),
            extra_packages=tuple(str(item) for item in payload.get("extra_packages", []) or []),
        )
        for name, payload in (variants_payload.get("variants", {}) or {}).items()
        if isinstance(payload, dict)
    }

    return DockerMetadata(
        root=root,
        distributions=distributions,
        package_groups=package_groups,
        layers=layers,
        variants=variants,
    )


def resolve_distribution_version(
    metadata: DockerMetadata,
    *,
    family: str,
    name: str,
) -> DockerDistributionVersion:
    key = f"{family}:{name}"
    try:
        return metadata.distributions[key]
    except KeyError as exc:
        known = ", ".join(sorted(metadata.distributions))
        raise ValueError(f"Unknown Docker distribution version: {key}. Known targets: {known}") from exc


def packages_for_variant(metadata: DockerMetadata, *, variant_name: str) -> list[str]:
    try:
        variant = metadata.variants[variant_name]
    except KeyError as exc:
        known = ", ".join(sorted(metadata.variants))
        raise ValueError(f"Unknown Docker variant: {variant_name}. Known variants: {known}") from exc

    packages: list[str] = []
    for layer_name in variant.layers:
        try:
            group_names = metadata.layers[layer_name]
        except KeyError as exc:
            known = ", ".join(sorted(metadata.layers))
            raise ValueError(f"Unknown Docker layer: {layer_name}. Known layers: {known}") from exc
        for group_name in group_names:
            try:
                group_packages = metadata.package_groups[group_name]
            except KeyError as exc:
                known = ", ".join(sorted(metadata.package_groups))
                raise ValueError(
                    f"Unknown Docker package group: {group_name}. Known groups: {known}"
                ) from exc
            packages.extend(group_packages)
    packages.extend(variant.extra_packages)

    seen: set[str] = set()
    ordered: list[str] = []
    for package in packages:
        if package in seen:
            continue
        seen.add(package)
        ordered.append(package)
    return ordered
