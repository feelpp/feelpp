from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


SEMVER_RE = re.compile(
    r"^(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)"
    r"(?:-(?P<prerelease>[0-9A-Za-z.-]+))?"
    r"(?:\+(?P<build>[0-9A-Za-z.-]+))?$"
)


@dataclass(frozen=True)
class SemanticVersion:
    major: int
    minor: int
    patch: int
    prerelease: str | None = None
    build: str | None = None

    @classmethod
    def parse(cls, raw: str) -> "SemanticVersion":
        match = SEMVER_RE.match(raw.strip())
        if not match:
            raise ValueError(f"Invalid semantic version: {raw}")
        return cls(
            major=int(match.group("major")),
            minor=int(match.group("minor")),
            patch=int(match.group("patch")),
            prerelease=match.group("prerelease") or None,
            build=match.group("build") or None,
        )

    @classmethod
    def from_debian_upstream(cls, raw: str) -> "SemanticVersion":
        semver_text = raw
        if "~" in semver_text:
            base, prerelease = semver_text.split("~", 1)
            semver_text = f"{base}-{prerelease}"
        return cls.parse(semver_text)

    @property
    def base(self) -> str:
        return f"{self.major}.{self.minor}.{self.patch}"

    @property
    def cmake_prerelease(self) -> str:
        if not self.prerelease:
            return ""
        return f"-{self.prerelease}"

    @property
    def debian_upstream(self) -> str:
        value = self.base
        if self.prerelease:
            value = f"{value}~{self.prerelease}"
        if self.build:
            value = f"{value}+{self.build}"
        return value

    def package_version(self, *, revision: int | str = 1, epoch: str | None = None) -> str:
        prefix = f"{epoch}:" if epoch else ""
        return f"{prefix}{self.debian_upstream}-{revision}"

    def __str__(self) -> str:
        value = self.base
        if self.prerelease:
            value = f"{value}-{self.prerelease}"
        if self.build:
            value = f"{value}+{self.build}"
        return value


@dataclass(frozen=True)
class DebianPackageVersion:
    upstream: str
    revision: str
    epoch: str | None = None

    @classmethod
    def parse(cls, raw: str) -> "DebianPackageVersion":
        candidate = raw.strip()
        epoch = None
        remainder = candidate
        if ":" in remainder:
            possible_epoch, rest = remainder.split(":", 1)
            if possible_epoch.isdigit():
                epoch = possible_epoch
                remainder = rest
        if "-" not in remainder:
            raise ValueError(f"Expected Debian package version with revision: {raw}")
        upstream, revision = remainder.rsplit("-", 1)
        if not upstream or not revision:
            raise ValueError(f"Invalid Debian package version: {raw}")
        return cls(upstream=upstream, revision=revision, epoch=epoch)

    @property
    def semver(self) -> SemanticVersion:
        return SemanticVersion.from_debian_upstream(self.upstream)

    def with_upstream(self, upstream: str, *, revision: int | str | None = None) -> "DebianPackageVersion":
        return DebianPackageVersion(
            upstream=upstream,
            revision=str(self.revision if revision is None else revision),
            epoch=self.epoch,
        )

    def increment_revision(self) -> "DebianPackageVersion":
        if not self.revision.isdigit():
            raise ValueError(f"Cannot increment non-numeric Debian revision: {self.revision}")
        return DebianPackageVersion(
            upstream=self.upstream,
            revision=str(int(self.revision) + 1),
            epoch=self.epoch,
        )

    def tag_safe(self) -> str:
        return str(self).replace(":", "-")

    def __str__(self) -> str:
        prefix = f"{self.epoch}:" if self.epoch else ""
        return f"{prefix}{self.upstream}-{self.revision}"


@dataclass(frozen=True)
class CMakeVersionRecord:
    name: str
    path: Path
    version: SemanticVersion

    def as_dict(self) -> dict[str, str]:
        return {"name": self.name, "path": str(self.path), "version": str(self.version)}


@dataclass(frozen=True)
class DebianPackageRecord:
    component: str
    dist: str
    source_name: str
    path: Path
    version: DebianPackageVersion
    flavor: str | None = None
    distribution: str | None = None
    urgency: str | None = None
    origin: str = "manifest"

    def as_dict(self) -> dict[str, str]:
        return {
            "component": self.component,
            "dist": self.dist,
            "flavor": self.flavor,
            "source_name": self.source_name,
            "path": str(self.path),
            "version": str(self.version),
            "distribution": self.distribution,
            "urgency": self.urgency,
            "origin": self.origin,
        }


@dataclass(frozen=True)
class MaintainerIdentity:
    name: str
    email: str

    @property
    def formatted(self) -> str:
        return f"{self.name} <{self.email}>"


@dataclass(frozen=True)
class PackageAvailabilityCheck:
    component: str
    dist: str
    flavor: str
    package_name: str
    expected_version: str
    url: str

    def as_dict(self) -> dict[str, str]:
        return {
            "component": self.component,
            "dist": self.dist,
            "flavor": self.flavor,
            "package_name": self.package_name,
            "expected_version": self.expected_version,
            "url": self.url,
        }


@dataclass(frozen=True)
class ContainerAvailabilityCheck:
    dist: str
    artifact_type: str
    candidate_refs: tuple[str, ...]

    def as_dict(self) -> dict[str, str]:
        return {
            "dist": self.dist,
            "artifact_type": self.artifact_type,
            "candidate_refs": list(self.candidate_refs),
        }


@dataclass(frozen=True)
class RepoVersionState:
    repo_root: Path
    cmake_versions: tuple[CMakeVersionRecord, ...]
    package_versions: tuple[DebianPackageRecord, ...]
    changelog_versions: tuple[DebianPackageRecord, ...] = ()

    def canonical_upstream_version(self) -> SemanticVersion:
        unique = {str(record.version): record.version for record in self.cmake_versions}
        if len(unique) != 1:
            details = ", ".join(f"{record.name}={record.version}" for record in self.cmake_versions)
            raise ValueError(f"Mismatched CMake versions: {details}")
        return next(iter(unique.values()))

    def require_matching_package_upstreams(self, expected: SemanticVersion | None = None) -> SemanticVersion:
        unique = {str(record.version.semver): record.version.semver for record in self.package_versions}
        if len(unique) != 1:
            details = ", ".join(
                f"{record.component}/{record.dist}={record.version}" for record in self.package_versions
            )
            raise ValueError(f"Mismatched package upstream versions: {details}")
        resolved = next(iter(unique.values()))
        if expected and str(resolved) != str(expected):
            raise ValueError(
                f"Package upstream version {resolved} does not match CMake version {expected}"
            )
        return resolved

    def package_record(self, component: str, dist: str) -> DebianPackageRecord:
        for record in self.package_versions:
            if record.component == component and record.dist == dist:
                return record
        raise KeyError(f"Missing package version for {component}/{dist}")

    def changelog_record(self, component: str, dist: str) -> DebianPackageRecord:
        for record in self.changelog_versions:
            if record.component == component and record.dist == dist:
                return record
        raise KeyError(f"Missing changelog version for {component}/{dist}")

    def as_dict(self) -> dict[str, object]:
        cmake_versions = [record.as_dict() for record in self.cmake_versions]
        package_versions = [record.as_dict() for record in self.package_versions]
        changelog_versions = [record.as_dict() for record in self.changelog_versions]
        cmake_consistent = len({record["version"] for record in cmake_versions}) == 1
        package_upstreams = {
            str(DebianPackageVersion.parse(record["version"]).semver) for record in package_versions
        }
        changelog_sync = True
        changelog_map = {
            (record.component, record.dist): str(record.version)
            for record in self.changelog_versions
        }
        for record in self.package_versions:
            if changelog_map.get((record.component, record.dist)) != str(record.version):
                changelog_sync = False
                break
        return {
            "repo_root": str(self.repo_root),
            "cmake_versions": cmake_versions,
            "package_versions": package_versions,
            "changelog_versions": changelog_versions,
            "consistency": {
                "cmake_versions": cmake_consistent,
                "package_upstreams": len(package_upstreams) == 1,
                "changelog_sync": changelog_sync,
            },
        }


@dataclass(frozen=True)
class ReleasePlan:
    requested_version: str
    release_kind: str
    tag: str
    title: str
    previous_tag: str | None
    repo_slug: str
    branch: str
    channel: str
    head_sha: str
    prerelease: bool
    package_notes: str
    generated_notes_preview: str
    package_checks: tuple[PackageAvailabilityCheck, ...]
    container_checks: tuple[ContainerAvailabilityCheck, ...]
    dry_run: bool

    def as_dict(self) -> dict[str, object]:
        return {
            "requested_version": self.requested_version,
            "release_kind": self.release_kind,
            "tag": self.tag,
            "title": self.title,
            "previous_tag": self.previous_tag,
            "repo_slug": self.repo_slug,
            "branch": self.branch,
            "channel": self.channel,
            "head_sha": self.head_sha,
            "prerelease": self.prerelease,
            "package_notes": self.package_notes,
            "generated_notes_preview": self.generated_notes_preview,
            "package_checks": [check.as_dict() for check in self.package_checks],
            "container_checks": [check.as_dict() for check in self.container_checks],
            "dry_run": self.dry_run,
        }
