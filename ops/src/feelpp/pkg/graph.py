from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib


SKIP_TEXT_TO_COMPONENT = {
    "skip feelpp": "feelpp",
    "skip toolboxes": "feelpp-toolboxes",
    "skip mor": "feelpp-mor",
}


@dataclass(frozen=True)
class ComponentSpec:
    name: str
    distros: tuple[str, ...]
    dependencies: tuple[str, ...]
    python_packages: tuple[str, ...]
    publish: bool
    package_revision: str = "1"
    package_epoch: str | None = None
    package_revision_by_dist: dict[str, str] | None = None

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "distros": list(self.distros),
            "dependencies": list(self.dependencies),
            "python_packages": list(self.python_packages),
            "publish": self.publish,
            "package_revision": self.package_revision,
            "package_epoch": self.package_epoch,
            "package_revision_by_dist": dict(self.package_revision_by_dist or {}),
        }

    def revision_for_dist(self, dist: str) -> str:
        overrides = self.package_revision_by_dist or {}
        return str(overrides.get(dist, self.package_revision))

    def package_version(self, upstream: str, dist: str) -> str:
        prefix = f"{self.package_epoch}:" if self.package_epoch else ""
        return f"{prefix}{upstream}-{self.revision_for_dist(dist)}"


@dataclass(frozen=True)
class Manifest:
    path: Path
    version: int
    default_components: tuple[str, ...]
    components: dict[str, ComponentSpec]


@dataclass(frozen=True)
class SkipSelection:
    components: frozenset[str]
    publish: bool


@dataclass(frozen=True)
class BuildPlan:
    dist: str
    components: tuple[ComponentSpec, ...]
    skipped_components: tuple[str, ...]
    publish_enabled: bool
    skip_text: str

    def as_dict(self) -> dict[str, object]:
        return {
            "dist": self.dist,
            "components": [component.as_dict() for component in self.components],
            "skipped_components": list(self.skipped_components),
            "publish_enabled": self.publish_enabled,
            "skip_text": self.skip_text,
        }


def load_manifest(path: Path) -> Manifest:
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    components_data = data.get("components", {})
    components: dict[str, ComponentSpec] = {}
    for name, component in components_data.items():
        components[name] = ComponentSpec(
            name=name,
            distros=tuple(component.get("distros", ())),
            dependencies=tuple(component.get("dependencies", ())),
            python_packages=tuple(component.get("python_packages", ())),
            publish=bool(component.get("publish", True)),
            package_revision=str(component.get("package_revision", "1")),
            package_epoch=(
                str(component["package_epoch"])
                if component.get("package_epoch") is not None
                else None
            ),
            package_revision_by_dist=(
                {
                    str(dist): str(revision)
                    for dist, revision in component.get("package_revision_by_dist", {}).items()
                }
                if component.get("package_revision_by_dist")
                else {}
            ),
        )
    return Manifest(
        path=path,
        version=int(data.get("version", 1)),
        default_components=tuple(data.get("default_components", ())),
        components=components,
    )


def parse_skip_text(text: str | None) -> SkipSelection:
    normalized = (text or "").lower()
    components = {
        component
        for token, component in SKIP_TEXT_TO_COMPONENT.items()
        if token in normalized
    }
    return SkipSelection(components=frozenset(components), publish="skip publish" in normalized)


def build_plan(
    manifest: Manifest,
    *,
    dist: str,
    requested_components: list[str] | None = None,
    skipped_components: set[str] | None = None,
    skip_text: str | None = None,
    publish: bool = False,
) -> BuildPlan:
    skip_selection = parse_skip_text(skip_text)
    skipped = set(skipped_components or ())
    skipped.update(skip_selection.components)

    ordered_names = list(requested_components or manifest.default_components)
    components: list[ComponentSpec] = []
    for name in ordered_names:
        if name not in manifest.components:
            raise ValueError(f"Unknown component: {name}")
        if name in skipped:
            continue
        component = manifest.components[name]
        if component.distros and dist not in component.distros:
            raise ValueError(f"Component {name} is not enabled for dist {dist}")
        components.append(component)

    return BuildPlan(
        dist=dist,
        components=tuple(components),
        skipped_components=tuple(sorted(skipped)),
        publish_enabled=publish and not skip_selection.publish,
        skip_text=skip_text or "",
    )
