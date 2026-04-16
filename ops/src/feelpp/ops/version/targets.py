from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import json
import re

from .models import (
    CMakeVersionRecord,
    DebianPackageRecord,
    DebianPackageVersion,
    MaintainerIdentity,
    MetadataVersionRecord,
    SemanticVersion,
)


def _replace_assignment(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf'(set\({re.escape(name)}\s+")([^"]*)("\s*\))')
    updated, count = pattern.subn(rf"\g<1>{value}\g<3>", text, count=1)
    if count != 1:
        raise ValueError(f"Unable to update {name}")
    return updated


def _read_assignment(text: str, name: str) -> str:
    pattern = re.compile(rf'set\({re.escape(name)}\s+"([^"]*)"\s*\)')
    match = pattern.search(text)
    if not match:
        raise ValueError(f"Unable to read {name}")
    return match.group(1)


@dataclass(frozen=True)
class CMakeVersionTarget:
    name: str
    path: Path

    def read(self) -> CMakeVersionRecord:
        text = self.path.read_text(encoding="utf-8")
        prerelease = _read_assignment(text, "FEELPP_VERSION_PRERELEASE") or None
        if prerelease and prerelease.startswith("-"):
            prerelease = prerelease[1:]
        return CMakeVersionRecord(
            name=self.name,
            path=self.path,
            version=SemanticVersion(
                major=int(_read_assignment(text, "FEELPP_VERSION_MAJOR")),
                minor=int(_read_assignment(text, "FEELPP_VERSION_MINOR")),
                patch=int(_read_assignment(text, "FEELPP_VERSION_MICRO")),
                prerelease=prerelease or None,
            ),
        )

    def write(self, version: SemanticVersion) -> None:
        text = self.path.read_text(encoding="utf-8")
        text = _replace_assignment(text, "FEELPP_VERSION_MAJOR", str(version.major))
        text = _replace_assignment(text, "FEELPP_VERSION_MINOR", str(version.minor))
        text = _replace_assignment(text, "FEELPP_VERSION_MICRO", str(version.patch))
        text = _replace_assignment(text, "FEELPP_VERSION_PRERELEASE", version.cmake_prerelease)
        self.path.write_text(text, encoding="utf-8")


def _coerce_metadata_version(raw: object, *, source: Path) -> SemanticVersion:
    text = str(raw or "").strip().strip('"').strip("'")
    if not text:
        raise ValueError(f"Unable to read version metadata from {source}")
    if text.startswith("v"):
        text = text[1:]
    return SemanticVersion.parse(text)


@dataclass(frozen=True)
class JsonMetadataVersionTarget:
    name: str
    path: Path
    key: str = "version"

    def read(self) -> MetadataVersionRecord:
        payload = json.loads(self.path.read_text(encoding="utf-8"))
        return MetadataVersionRecord(
            name=self.name,
            path=self.path,
            version=_coerce_metadata_version(payload.get(self.key), source=self.path),
        )

    def write(self, version: SemanticVersion) -> None:
        payload = json.loads(self.path.read_text(encoding="utf-8"))
        payload[self.key] = f"v{version}"
        self.path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _replace_or_insert_line(text: str, *, key: str, value: str, insert_after: str | None = None) -> str:
    pattern = re.compile(rf"^(?P<prefix>{re.escape(key)}:\s*)(?P<value>.*)$", re.MULTILINE)
    replacement = rf"\g<prefix>{value}"
    updated, count = pattern.subn(replacement, text, count=1)
    if count == 1:
        return updated

    lines = text.splitlines()
    insert_at = len(lines)
    if insert_after is not None:
        for index, line in enumerate(lines):
            if line.startswith(f"{insert_after}:"):
                insert_at = index + 1
                break
    lines.insert(insert_at, f"{key}: {value}")
    rendered = "\n".join(lines)
    if text.endswith("\n") or not text:
        rendered += "\n"
    return rendered


@dataclass(frozen=True)
class CitationCffVersionTarget:
    name: str
    path: Path

    def read(self) -> MetadataVersionRecord:
        text = self.path.read_text(encoding="utf-8")
        match = re.search(r"^version:\s*(?P<value>.+?)\s*$", text, re.MULTILINE)
        if not match:
            raise ValueError(f"Unable to read version metadata from {self.path}")
        return MetadataVersionRecord(
            name=self.name,
            path=self.path,
            version=_coerce_metadata_version(match.group("value"), source=self.path),
        )

    def write(self, version: SemanticVersion) -> None:
        text = self.path.read_text(encoding="utf-8")
        updated = _replace_or_insert_line(text, key="version", value=f"v{version}", insert_after="title")
        self.path.write_text(updated, encoding="utf-8")


CHANGELOG_HEADER_RE = re.compile(
    r"^(?P<source>\S+) \((?P<version>[^)]+)\) (?P<distribution>\S+); urgency=(?P<urgency>\S+)\s*$"
)


@dataclass(frozen=True)
class DebianChangelogTarget:
    component: str
    dist: str
    path: Path

    def read(self) -> DebianPackageRecord:
        text = self.path.read_text(encoding="utf-8")
        first_line = text.splitlines()[0]
        match = CHANGELOG_HEADER_RE.match(first_line)
        if not match:
            raise ValueError(f"Unsupported changelog header in {self.path}")
        return DebianPackageRecord(
            component=self.component,
            dist=self.dist,
            source_name=match.group("source"),
            path=self.path,
            version=DebianPackageVersion.parse(match.group("version")),
            flavor=None,
            distribution=match.group("distribution"),
            urgency=match.group("urgency"),
            origin="changelog",
        )

    def prepend_entry(
        self,
        version: DebianPackageVersion,
        *,
        message: str,
        identity: MaintainerIdentity,
        timestamp: datetime,
    ) -> DebianPackageRecord:
        current = self.read()
        if str(current.version) == str(version):
            return current

        text = self.path.read_text(encoding="utf-8")
        rendered_timestamp = timestamp.strftime("%a, %d %b %Y %H:%M:%S %z")
        entry = (
            f"{current.source_name} ({version}) {current.distribution}; urgency={current.urgency}\n\n"
            f"  * {message}\n\n"
            f" -- {identity.formatted}  {rendered_timestamp}\n\n"
        )
        self.path.write_text(entry + text, encoding="utf-8")
        return DebianPackageRecord(
            component=current.component,
            dist=current.dist,
            source_name=current.source_name,
            path=current.path,
            version=version,
            flavor=current.flavor,
            distribution=current.distribution,
            urgency=current.urgency,
            origin=current.origin,
        )
