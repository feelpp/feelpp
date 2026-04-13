from __future__ import annotations

from datetime import datetime
from pathlib import Path
import os
from typing import Callable, TypeVar

from feelpp.pkg.config import detect_flavor, discover_repo_root
from feelpp.pkg.graph import ComponentSpec, load_manifest
from feelpp.pkg.shell import run_capture

from .models import DebianPackageRecord, DebianPackageVersion, MaintainerIdentity, RepoVersionState, SemanticVersion
from .targets import CMakeVersionTarget, DebianChangelogTarget


CMAKE_TARGET_PATHS = (
    ("feelpp", Path("feelpp.version.cmake")),
)

T = TypeVar("T")


def current_timestamp() -> datetime:
    return datetime.now().astimezone()


def resolve_maintainer_identity(repo_root: Path) -> MaintainerIdentity:
    name = os.getenv("DEBFULLNAME") or _git_config(repo_root, "user.name") or "Feel++ Release Automation"
    email = (
        os.getenv("DEBEMAIL")
        or os.getenv("EMAIL")
        or _git_config(repo_root, "user.email")
        or "noreply@feelpp.org"
    )
    return MaintainerIdentity(name=name, email=email)


def _git_config(repo_root: Path, key: str) -> str | None:
    try:
        value = run_capture(["git", "config", "--get", key], cwd=repo_root, check=False).strip()
    except FileNotFoundError:
        return None
    return value or None


class VersionRepository:
    def __init__(self, repo_root: str | Path | None = None) -> None:
        self.repo_root = Path(repo_root).expanduser().resolve() if repo_root else discover_repo_root()

    @property
    def manifest_path(self) -> Path:
        return self.repo_root / "packaging" / "manifest" / "components.toml"

    def cmake_targets(self) -> tuple[CMakeVersionTarget, ...]:
        return tuple(
            CMakeVersionTarget(name=name, path=self.repo_root / relative_path)
            for name, relative_path in CMAKE_TARGET_PATHS
        )

    def changelog_targets(self) -> tuple[DebianChangelogTarget, ...]:
        base = self.repo_root / "packaging" / "debian"
        targets: list[DebianChangelogTarget] = []
        for path in sorted(base.glob("*/*/debian/changelog")):
            component = path.parts[-4]
            dist = path.parts[-3]
            targets.append(DebianChangelogTarget(component=component, dist=dist, path=path))
        return tuple(targets)

    def manifest(self):
        return load_manifest(self.manifest_path)

    def package_versions(self, version: SemanticVersion) -> tuple[DebianPackageRecord, ...]:
        manifest = self.manifest()
        records: list[DebianPackageRecord] = []
        for component_name in sorted(manifest.components):
            component = manifest.components[component_name]
            for dist in component.distros:
                records.append(
                    DebianPackageRecord(
                        component=component.name,
                        dist=dist,
                        flavor=detect_flavor(dist),
                        source_name=component.name,
                        path=self.manifest_path,
                        version=DebianPackageVersion.parse(component.package_version(version.debian_upstream, dist)),
                        origin="manifest",
                    )
                )
        return tuple(records)

    def _mutable_paths(self) -> tuple[Path, ...]:
        unique: dict[Path, None] = {}
        for target in self.cmake_targets():
            unique[target.path] = None
        unique[self.manifest_path] = None
        for target in self.changelog_targets():
            unique[target.path] = None
        return tuple(unique.keys())

    def _run_with_rollback(self, operation: Callable[[], T]) -> T:
        snapshots = {
            path: path.read_text(encoding="utf-8")
            for path in self._mutable_paths()
            if path.exists()
        }
        try:
            return operation()
        finally:
            for path, content in snapshots.items():
                path.write_text(content, encoding="utf-8")

    def _component_section_headers(self, component: str) -> tuple[str, ...]:
        return (
            f"[components.{component}]",
            f'[components."{component}"]',
        )

    def _render_manifest_value(self, value: str | dict[str, str]) -> str:
        if isinstance(value, dict):
            if not value:
                return "{}"
            rendered = ", ".join(
                f'{dist} = "{revision}"' for dist, revision in sorted(value.items())
            )
            return f"{{ {rendered} }}"
        return f'"{value}"'

    def _set_component_manifest_field(
        self,
        component: str,
        field: str,
        value: str | dict[str, str] | None,
    ) -> None:
        path = self.manifest_path
        lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
        headers = self._component_section_headers(component)

        start_index: int | None = None
        for index, raw_line in enumerate(lines):
            if raw_line.strip() in headers:
                start_index = index
                break
        if start_index is None:
            raise KeyError(f"Component {component} not found in {path}")

        end_index = len(lines)
        for index in range(start_index + 1, len(lines)):
            stripped = lines[index].strip()
            if stripped.startswith("[components.") or (
                stripped.startswith("[") and not stripped.startswith("[components.")
            ):
                end_index = index
                break

        for index in range(start_index + 1, end_index):
            if lines[index].lstrip().startswith(f"{field} ="):
                if value is None:
                    del lines[index]
                else:
                    lines[index] = f"{field} = {self._render_manifest_value(value)}\n"
                path.write_text("".join(lines), encoding="utf-8")
                return

        if value is None:
            return

        insert_at = end_index
        while insert_at > start_index + 1 and lines[insert_at - 1].strip() == "":
            insert_at -= 1
        lines.insert(insert_at, f"{field} = {self._render_manifest_value(value)}\n")
        path.write_text("".join(lines), encoding="utf-8")

    def _set_component_package_revision(self, component: str, revision: str) -> None:
        self._set_component_manifest_field(component, "package_revision", revision)

    def _set_component_package_revision_by_dist(
        self,
        component: str,
        revisions: dict[str, str] | None,
    ) -> None:
        self._set_component_manifest_field(
            component,
            "package_revision_by_dist",
            revisions or None,
        )

    def _selected_changelog_targets(
        self,
        dists: tuple[str, ...] | None = None,
    ) -> tuple[DebianChangelogTarget, ...]:
        if not dists:
            return self.changelog_targets()
        selected = set(dists)
        return tuple(target for target in self.changelog_targets() if target.dist in selected)

    def _normalize_component_revisions(
        self,
        component: ComponentSpec,
        per_dist_revisions: dict[str, str],
    ) -> tuple[str, dict[str, str]]:
        current_default = str(component.package_revision)
        enabled_dists = [dist for dist in component.distros if dist in per_dist_revisions]
        if not enabled_dists:
            return current_default, {}

        distinct = {per_dist_revisions[dist] for dist in enabled_dists}
        if len(distinct) == 1:
            return next(iter(distinct)), {}

        counts: dict[str, int] = {}
        for dist in enabled_dists:
            revision = per_dist_revisions[dist]
            counts[revision] = counts.get(revision, 0) + 1
        max_count = max(counts.values())
        candidates = {revision for revision, count in counts.items() if count == max_count}
        default_revision = current_default if current_default in candidates else sorted(candidates)[0]
        overrides = {
            dist: revision
            for dist, revision in per_dist_revisions.items()
            if revision != default_revision
        }
        return default_revision, overrides

    def _apply_component_dist_revision_bump(
        self,
        component: ComponentSpec,
        *,
        upstream: SemanticVersion,
        dists: tuple[str, ...] | None = None,
    ) -> None:
        selected_dists = set(dists or component.distros)
        selected_dists &= set(component.distros)
        if not selected_dists:
            return

        per_dist_revisions = {
            dist: component.revision_for_dist(dist)
            for dist in component.distros
        }
        for dist in selected_dists:
            current = DebianPackageVersion.parse(component.package_version(upstream.debian_upstream, dist))
            per_dist_revisions[dist] = current.increment_revision().revision

        default_revision, overrides = self._normalize_component_revisions(component, per_dist_revisions)
        self._set_component_package_revision(component.name, default_revision)
        self._set_component_package_revision_by_dist(component.name, overrides)

    def _reset_component_package_revisions(self, component: str) -> None:
        self._set_component_package_revision(component, "1")
        self._set_component_package_revision_by_dist(component, None)

    def sync_changelogs(
        self,
        *,
        message: str = "Packaging metadata sync",
        identity: MaintainerIdentity | None = None,
        timestamp: datetime | None = None,
        dists: tuple[str, ...] | None = None,
        dry_run: bool = False,
    ) -> RepoVersionState:
        def apply() -> RepoVersionState:
            state = self.read_state()
            canonical = state.canonical_upstream_version()
            manifest = self.manifest()
            resolved_identity = identity or resolve_maintainer_identity(self.repo_root)
            resolved_timestamp = timestamp or current_timestamp()

            for target in self._selected_changelog_targets(dists):
                component = manifest.components[target.component]
                target.prepend_entry(
                    DebianPackageVersion.parse(component.package_version(canonical.debian_upstream, target.dist)),
                    message=message,
                    identity=resolved_identity,
                    timestamp=resolved_timestamp,
                )
            return self.read_state()

        return self._run_with_rollback(apply) if dry_run else apply()

    def read_state(self) -> RepoVersionState:
        cmake_versions = tuple(target.read() for target in self.cmake_targets())
        canonical = cmake_versions[0].version
        return RepoVersionState(
            repo_root=self.repo_root,
            cmake_versions=cmake_versions,
            package_versions=self.package_versions(canonical),
            changelog_versions=tuple(target.read() for target in self.changelog_targets()),
        )

    def bump_upstream(
        self,
        version: SemanticVersion,
        *,
        message: str = "New upstream release",
        identity: MaintainerIdentity | None = None,
        timestamp: datetime | None = None,
        dry_run: bool = False,
    ) -> RepoVersionState:
        def apply() -> RepoVersionState:
            resolved_identity = identity or resolve_maintainer_identity(self.repo_root)
            resolved_timestamp = timestamp or current_timestamp()
            manifest = self.manifest()

            for target in self.cmake_targets():
                target.write(version)

            for component_name in manifest.components:
                self._reset_component_package_revisions(component_name)

            refreshed_manifest = self.manifest()

            for target in self.changelog_targets():
                target.prepend_entry(
                    DebianPackageVersion.parse(
                        refreshed_manifest.components[target.component].package_version(
                            version.debian_upstream,
                            target.dist,
                        )
                    ),
                    message=message,
                    identity=resolved_identity,
                    timestamp=resolved_timestamp,
                )
            return self.read_state()

        return self._run_with_rollback(apply) if dry_run else apply()

    def bump_revision(
        self,
        *,
        message: str = "Packaging revision update",
        identity: MaintainerIdentity | None = None,
        timestamp: datetime | None = None,
        dists: tuple[str, ...] | None = None,
        dry_run: bool = False,
    ) -> RepoVersionState:
        def apply() -> RepoVersionState:
            state = self.read_state()
            canonical = state.canonical_upstream_version()
            state.require_matching_package_upstreams(expected=canonical)
            manifest = self.manifest()

            resolved_identity = identity or resolve_maintainer_identity(self.repo_root)
            resolved_timestamp = timestamp or current_timestamp()

            for component in manifest.components.values():
                self._apply_component_dist_revision_bump(
                    component,
                    upstream=canonical,
                    dists=dists,
                )

            return self.sync_changelogs(
                message=message,
                identity=resolved_identity,
                timestamp=resolved_timestamp,
                dists=dists,
            )

        return self._run_with_rollback(apply) if dry_run else apply()
