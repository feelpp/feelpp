from __future__ import annotations

from datetime import datetime
from pathlib import Path
import os

from feelpp.pkg.config import discover_repo_root
from feelpp.pkg.shell import run_capture

from .models import MaintainerIdentity, RepoVersionState, SemanticVersion
from .targets import CMakeVersionTarget, DebianChangelogTarget


CMAKE_TARGET_PATHS = (
    ("feelpp", Path("feelpp.version.cmake")),
    ("feelpp-toolboxes", Path("toolboxes") / "cmake" / "feelpp.version.cmake"),
    ("feelpp-mor", Path("mor") / "cmake" / "feelpp.version.cmake"),
)


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

    def read_state(self) -> RepoVersionState:
        return RepoVersionState(
            repo_root=self.repo_root,
            cmake_versions=tuple(target.read() for target in self.cmake_targets()),
            package_versions=tuple(target.read() for target in self.changelog_targets()),
        )

    def bump_upstream(
        self,
        version: SemanticVersion,
        *,
        message: str = "New upstream release",
        identity: MaintainerIdentity | None = None,
        timestamp: datetime | None = None,
    ) -> RepoVersionState:
        resolved_identity = identity or resolve_maintainer_identity(self.repo_root)
        resolved_timestamp = timestamp or current_timestamp()

        for target in self.cmake_targets():
            target.write(version)

        for target in self.changelog_targets():
            current = target.read()
            target.prepend_entry(
                current.version.with_upstream(version.debian_upstream, revision=1),
                message=message,
                identity=resolved_identity,
                timestamp=resolved_timestamp,
            )
        return self.read_state()

    def bump_revision(
        self,
        *,
        message: str = "Packaging revision update",
        identity: MaintainerIdentity | None = None,
        timestamp: datetime | None = None,
    ) -> RepoVersionState:
        state = self.read_state()
        canonical = state.canonical_upstream_version()
        state.require_matching_package_upstreams(expected=canonical)

        resolved_identity = identity or resolve_maintainer_identity(self.repo_root)
        resolved_timestamp = timestamp or current_timestamp()

        for target in self.changelog_targets():
            current = target.read()
            target.prepend_entry(
                current.version.increment_revision(),
                message=message,
                identity=resolved_identity,
                timestamp=resolved_timestamp,
            )
        return self.read_state()
