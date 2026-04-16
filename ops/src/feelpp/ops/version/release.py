from __future__ import annotations

from dataclasses import replace
from functools import cmp_to_key
from pathlib import Path
import gzip
import json
import lzma
import os
import re
import subprocess
import tempfile
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from feelpp.pkg.config import detect_channel, detect_flavor
from feelpp.pkg.graph import load_manifest
from feelpp.pkg.image import normalize_package_version_tag
from feelpp.pkg.shell import run_capture, run_checked

from .publications import HalPublicationService
from .models import (
    ContainerAvailabilityCheck,
    DebianPackageVersion,
    GitHubContributor,
    PackageAvailabilityCheck,
    ReleasePlan,
    SemanticVersion,
)
from .repository import VersionRepository


OCI_REPOSITORY = os.getenv("FEELPP_PKG_OCI_REPOSITORY") or "feelpp/feelpp"
OCI_REGISTRY = os.getenv("FEELPP_PKG_OCI_REGISTRY") or "ghcr.io"
APPTAINER_TAG_SUFFIXES = ("-sif", "_sif")
GENERATED_NOTES_UNAVAILABLE = "* GitHub-generated release notes preview unavailable."
GITHUB_LOGIN_RE = re.compile(r"(?<![A-Za-z0-9/])@(?P<login>[A-Za-z0-9][A-Za-z0-9-]*(?:\[[A-Za-z0-9-]+\])?)")


class ReleaseService:
    def __init__(self, repo_root: str | Path | None = None) -> None:
        self.repository = VersionRepository(repo_root=repo_root)
        self.repo_root = self.repository.repo_root

    def prepare_release(
        self,
        requested_version: str,
        *,
        dry_run: bool = False,
        dists: tuple[str, ...] | None = None,
        publication_rows: int | None = None,
        publication_since: str | None = None,
    ) -> ReleasePlan:
        state = self.repository.read_state()
        canonical_upstream = state.canonical_upstream_version()
        state.require_matching_package_upstreams(expected=canonical_upstream)
        selected_records = self._selected_package_records(state, dists)
        dist_versions = self._package_versions_by_dist(selected_records)
        selected_dists = sorted(dist_versions)
        omitted_dists = sorted({record.dist for record in state.package_versions} - set(selected_dists))

        try:
            requested_semver = SemanticVersion.parse(requested_version)
        except ValueError:
            requested_semver = None

        requested_package: DebianPackageVersion | None = None
        release_kind: str
        tag: str
        title: str
        prerelease: bool
        if requested_semver is not None:
            release_kind = "upstream"
            if str(canonical_upstream) != str(requested_semver):
                raise ValueError(
                    f"Requested version {requested_semver} does not match repository version {canonical_upstream}"
                )
            tag = f"v{requested_semver}"
            title = tag
            prerelease = requested_semver.prerelease is not None
        else:
            requested_package = DebianPackageVersion.parse(requested_version)
            release_kind = "packaging"
            if requested_package.semver != canonical_upstream:
                raise ValueError(
                    f"Requested package version {requested_package} does not match repository version {canonical_upstream}"
                )
            mismatched = [
                str(record.version) for record in selected_records if str(record.version) != str(requested_package)
            ]
            if mismatched:
                raise ValueError(
                    "Packaging release requires selected package versions to match the requested package version; "
                    f"found: {', '.join(mismatched)}"
                )
            tag = f"pkg/{requested_package.tag_safe()}"
            title = f"Packaging update {requested_package}"
            prerelease = True

        branch = self._current_branch()
        channel = detect_channel(branch)
        head_sha = self._git_capture(["rev-parse", "HEAD"]).strip()
        repo_slug = self._repo_slug()

        self._ensure_clean_worktree()
        self._ensure_head_pushed()
        self._ensure_tag_absent(tag)
        self._ensure_github_checks_green(repo_slug, head_sha)

        package_checks = self._package_checks(selected_records, channel=channel)
        for check in package_checks:
            self._ensure_package_available(check)

        container_checks = self._container_checks(dist_versions)
        for check in container_checks:
            self._ensure_container_available(check)

        previous_tag = self._previous_tag(tag)
        package_notes = self._package_notes(
            package_versions_by_dist=dist_versions,
            package_checks=package_checks,
            channel=channel,
            omitted_dists=omitted_dists,
        )
        publication_notes = self._publication_notes(
            rows=publication_rows,
            since=publication_since,
        )
        if publication_notes:
            package_notes = f"{package_notes}\n\n{publication_notes}"
        generated_notes_preview = self._generated_notes_preview(
            repo_slug=repo_slug,
            tag=tag,
            head_sha=head_sha,
            previous_tag=previous_tag,
        )
        return ReleasePlan(
            requested_version=requested_version,
            release_kind=release_kind,
            tag=tag,
            title=title,
            previous_tag=previous_tag,
            repo_slug=repo_slug,
            branch=branch,
            channel=channel,
            head_sha=head_sha,
            prerelease=prerelease,
            package_notes=package_notes,
            generated_notes_preview=generated_notes_preview,
            package_checks=tuple(package_checks),
            container_checks=tuple(container_checks),
            dry_run=dry_run,
        )

    def execute_release(
        self,
        requested_version: str,
        *,
        dry_run: bool = False,
        dists: tuple[str, ...] | None = None,
        publication_rows: int | None = None,
        publication_since: str | None = None,
    ) -> ReleasePlan:
        plan = self.prepare_release(
            requested_version,
            dry_run=dry_run,
            dists=dists,
            publication_rows=publication_rows,
            publication_since=publication_since,
        )
        if dry_run:
            return plan

        release_head_sha = self._sync_release_metadata(plan)
        run_checked(
            ["git", "tag", "-a", plan.tag, "-m", plan.title, release_head_sha],
            cwd=self.repo_root,
        )
        run_checked(["git", "push", "origin", f"refs/tags/{plan.tag}"], cwd=self.repo_root)

        with tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False) as handle:
            handle.write(self._render_release_notes(plan))
            notes_path = Path(handle.name)

        try:
            command = [
                "gh",
                "release",
                "create",
                plan.tag,
                "--verify-tag",
                "--title",
                plan.title,
                "--notes-file",
                str(notes_path),
            ]
            if plan.prerelease:
                command.append("--prerelease")
            run_checked(command, cwd=self.repo_root)
        finally:
            notes_path.unlink(missing_ok=True)
        return replace(plan, head_sha=release_head_sha)

    def _current_branch(self) -> str:
        return self._git_capture(["branch", "--show-current"]).strip() or "develop"

    def _repo_slug(self) -> str:
        remote = self._git_capture(["remote", "get-url", "origin"]).strip()
        if remote.endswith(".git"):
            remote = remote[:-4]
        if remote.startswith("git@github.com:"):
            return remote.split(":", 1)[1]
        if remote.startswith("https://github.com/"):
            return remote.split("https://github.com/", 1)[1]
        if remote.startswith("ssh://git@github.com/"):
            return remote.split("ssh://git@github.com/", 1)[1]
        raise ValueError(f"Unsupported GitHub remote URL: {remote}")

    def _git_capture(self, args: list[str], *, check: bool = True) -> str:
        return run_capture(["git", *args], cwd=self.repo_root, check=check)

    def _ensure_clean_worktree(self) -> None:
        status = self._git_capture(["status", "--short"]).strip()
        if status:
            raise RuntimeError("Release requires a clean git worktree")

    def _ensure_head_pushed(self) -> None:
        upstream = self._git_capture(["rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"]).strip()
        if not upstream:
            raise RuntimeError("Release requires the current branch to track a remote upstream")
        local_head = self._git_capture(["rev-parse", "HEAD"]).strip()
        remote_head = self._git_capture(["rev-parse", upstream]).strip()
        if local_head != remote_head:
            raise RuntimeError(f"Local HEAD {local_head} is not pushed to {upstream}")

    def _ensure_tag_absent(self, tag: str) -> None:
        existing = self._git_capture(["tag", "--list", tag], check=False).strip()
        if existing:
            raise RuntimeError(f"Tag already exists: {tag}")

    def _previous_tag(self, tag: str) -> str | None:
        if tag.startswith("pkg/"):
            previous = self._git_capture(
                ["describe", "--tags", "--abbrev=0", "--match", "pkg/*"],
                check=False,
            ).strip()
            if previous == tag:
                previous = self._git_capture(
                    ["describe", "--tags", "--abbrev=0", "--match", "pkg/*", f"{tag}^"],
                    check=False,
                ).strip()
            return previous or None

        if not tag.startswith("v"):
            return None

        current = SemanticVersion.parse(tag[1:])
        raw_tags = self._git_capture(["tag", "--merged", "HEAD", "--list", "v*"], check=False)
        candidates: list[tuple[str, SemanticVersion]] = []
        for raw_tag in raw_tags.splitlines():
            candidate_tag = raw_tag.strip()
            if not candidate_tag or candidate_tag == tag or not candidate_tag.startswith("v"):
                continue
            try:
                candidate_version = SemanticVersion.parse(candidate_tag[1:])
            except ValueError:
                continue
            candidates.append((candidate_tag, candidate_version))

        if not candidates:
            return None

        if current.prerelease is not None:
            same_series = [
                (candidate_tag, candidate_version)
                for candidate_tag, candidate_version in candidates
                if candidate_version.base == current.base
                and candidate_version.prerelease is not None
                and self._compare_semver(candidate_version, current) < 0
            ]
            if same_series:
                return max(same_series, key=cmp_to_key(self._compare_semver_tag))[0]

        stable_candidates = [
            (candidate_tag, candidate_version)
            for candidate_tag, candidate_version in candidates
            if candidate_version.prerelease is None
            and self._compare_semver(candidate_version, current) < 0
        ]
        if stable_candidates:
            return max(stable_candidates, key=cmp_to_key(self._compare_semver_tag))[0]
        return None

    def _compare_semver_tag(
        self,
        left: tuple[str, SemanticVersion],
        right: tuple[str, SemanticVersion],
    ) -> int:
        return self._compare_semver(left[1], right[1])

    def _compare_semver(self, left: SemanticVersion, right: SemanticVersion) -> int:
        left_core = (left.major, left.minor, left.patch)
        right_core = (right.major, right.minor, right.patch)
        if left_core < right_core:
            return -1
        if left_core > right_core:
            return 1
        return self._compare_prerelease(left.prerelease, right.prerelease)

    def _compare_prerelease(self, left: str | None, right: str | None) -> int:
        if left is None and right is None:
            return 0
        if left is None:
            return 1
        if right is None:
            return -1

        left_parts = left.split(".")
        right_parts = right.split(".")
        for left_part, right_part in zip(left_parts, right_parts):
            left_numeric = left_part.isdigit()
            right_numeric = right_part.isdigit()
            if left_numeric and right_numeric:
                left_value = int(left_part)
                right_value = int(right_part)
                if left_value < right_value:
                    return -1
                if left_value > right_value:
                    return 1
                continue
            if left_numeric and not right_numeric:
                return -1
            if not left_numeric and right_numeric:
                return 1
            if left_part < right_part:
                return -1
            if left_part > right_part:
                return 1

        if len(left_parts) < len(right_parts):
            return -1
        if len(left_parts) > len(right_parts):
            return 1
        return 0

    def _generated_notes_preview(
        self,
        *,
        repo_slug: str,
        tag: str,
        head_sha: str,
        previous_tag: str | None,
    ) -> str:
        command = [
            "gh",
            "api",
            f"repos/{repo_slug}/releases/generate-notes",
            "-X",
            "POST",
            "-f",
            f"tag_name={tag}",
            "-f",
            f"target_commitish={head_sha}",
        ]
        if previous_tag:
            command.extend(["-f", f"previous_tag_name={previous_tag}"])
        try:
            payload = json.loads(run_capture(command, cwd=self.repo_root))
        except (FileNotFoundError, RuntimeError, ValueError, json.JSONDecodeError):
            return GENERATED_NOTES_UNAVAILABLE
        body = str(payload.get("body", "")).strip()
        return body or GENERATED_NOTES_UNAVAILABLE

    def _release_semver(self, requested_version: str) -> SemanticVersion:
        try:
            return SemanticVersion.parse(requested_version)
        except ValueError:
            return DebianPackageVersion.parse(requested_version).semver

    def _metadata_path_args(self) -> list[str]:
        return [str(path.relative_to(self.repo_root)) for path in self.repository.metadata_paths()]

    def _metadata_paths_changed(self) -> bool:
        metadata_paths = self._metadata_path_args()
        if not metadata_paths:
            return False
        status = self._git_capture(["status", "--short", "--", *metadata_paths]).strip()
        return bool(status)

    def _sync_release_metadata(self, plan: ReleasePlan) -> str:
        release_version = self._release_semver(plan.requested_version)
        self.repository.sync_metadata(version=release_version, contributors=self._release_contributors(plan))
        if not self._metadata_paths_changed():
            return self._git_capture(["rev-parse", "HEAD"]).strip()

        metadata_paths = self._metadata_path_args()
        run_checked(["git", "add", "--", *metadata_paths], cwd=self.repo_root)
        run_checked(
            ["git", "commit", "-m", f"chore(release): sync metadata for {plan.tag} [ci skip]"],
            cwd=self.repo_root,
        )
        run_checked(["git", "push", "origin", plan.branch], cwd=self.repo_root)
        return self._git_capture(["rev-parse", "HEAD"]).strip()

    def _release_contributors(self, plan: ReleasePlan) -> tuple[GitHubContributor, ...]:
        if not plan.generated_notes_preview or plan.generated_notes_preview == GENERATED_NOTES_UNAVAILABLE:
            return ()

        logins: list[str] = []
        seen_logins: set[str] = set()
        for match in GITHUB_LOGIN_RE.finditer(plan.generated_notes_preview):
            login = str(match.group("login") or "").strip()
            if not login:
                continue
            normalized = login.lower()
            if normalized in seen_logins:
                continue
            seen_logins.add(normalized)
            logins.append(login)

        contributors: list[GitHubContributor] = []
        for login in logins:
            contributor = self._github_contributor(login)
            if contributor is not None:
                contributors.append(contributor)
        return tuple(contributors)

    def _github_contributor(self, login: str) -> GitHubContributor | None:
        try:
            payload = json.loads(run_capture(["gh", "api", f"users/{login}"], cwd=self.repo_root))
        except (FileNotFoundError, RuntimeError, ValueError, json.JSONDecodeError):
            return GitHubContributor(login=login, name=None)

        if str(payload.get("type") or "").strip().lower() == "bot":
            return None

        resolved_login = str(payload.get("login") or login).strip() or login
        resolved_name = str(payload.get("name") or "").strip() or None
        return GitHubContributor(login=resolved_login, name=resolved_name)

    def _render_release_notes(self, plan: ReleasePlan) -> str:
        sections = [plan.package_notes.strip()]
        if plan.generated_notes_preview and plan.generated_notes_preview != GENERATED_NOTES_UNAVAILABLE:
            sections.append(plan.generated_notes_preview.strip())
        return "\n\n".join(section for section in sections if section)

    def _ensure_github_checks_green(self, repo_slug: str, sha: str) -> None:
        status_output = run_capture(
            ["gh", "api", f"repos/{repo_slug}/commits/{sha}/status"],
            cwd=self.repo_root,
        )
        status_data = json.loads(status_output)
        statuses = status_data.get("statuses", [])
        if statuses and status_data.get("state") != "success":
            raise RuntimeError(f"GitHub commit status is not green for {sha}: {status_data.get('state')}")

        checks_output = run_capture(
            ["gh", "api", f"repos/{repo_slug}/commits/{sha}/check-runs"],
            cwd=self.repo_root,
        )
        checks_data = json.loads(checks_output)
        failures = [
            check["name"]
            for check in checks_data.get("check_runs", [])
            if check.get("status") != "completed"
            or check.get("conclusion") not in {"success", "neutral", "skipped"}
        ]
        if failures:
            raise RuntimeError(f"GitHub checks are not green for {sha}: {', '.join(failures)}")

    def _selected_package_records(self, state, dists: tuple[str, ...] | None) -> tuple:
        if not dists:
            return state.package_versions
        selected_dists = tuple(dict.fromkeys(dists))
        available = {record.dist for record in state.package_versions}
        missing = [dist for dist in selected_dists if dist not in available]
        if missing:
            raise ValueError(f"Unknown or unpublished dist selection: {', '.join(missing)}")
        return tuple(record for record in state.package_versions if record.dist in set(selected_dists))

    def _package_versions_by_dist(self, package_records) -> dict[str, DebianPackageVersion]:
        versions_by_dist: dict[str, DebianPackageVersion] = {}
        for record in package_records:
            existing = versions_by_dist.get(record.dist)
            if existing is None:
                versions_by_dist[record.dist] = record.version
                continue
            if str(existing) != str(record.version):
                raise RuntimeError(
                    f"Release requires a single package version per dist, but {record.dist} "
                    f"has both {existing} and {record.version}"
                )
        return versions_by_dist

    def _package_checks(self, package_records, *, channel: str) -> list[PackageAvailabilityCheck]:
        manifest = load_manifest(self.repository.manifest_path)
        checks: list[PackageAvailabilityCheck] = []
        for record in package_records:
            component = manifest.components[record.component]
            package_name = component.python_packages[0] if component.python_packages else record.component
            url = (
                f"http://apt.feelpp.org/{record.flavor}/{record.dist}/dists/{record.dist}/{channel}/binary-amd64"
            )
            checks.append(
                PackageAvailabilityCheck(
                    component=record.component,
                    dist=record.dist,
                    flavor=record.flavor,
                    package_name=package_name,
                    expected_version=str(record.version),
                    url=url,
                )
            )
        return checks

    def _ensure_package_available(self, check: PackageAvailabilityCheck) -> None:
        index_text = self._fetch_package_index(check.url)
        package_marker = f"Package: {check.package_name}\n"
        version_marker = f"Version: {check.expected_version}\n"
        if package_marker not in index_text or version_marker not in index_text:
            raise RuntimeError(
                f"Package {check.package_name} version {check.expected_version} was not found at {check.url}"
            )

    def _fetch_package_index(self, base_url: str) -> str:
        candidates = (
            f"{base_url}/Packages.xz",
            f"{base_url}/Packages.gz",
            f"{base_url}/Packages",
        )
        for candidate in candidates:
            try:
                with urlopen(candidate, timeout=10) as response:
                    payload = response.read()
            except (HTTPError, URLError):
                continue
            if candidate.endswith(".xz"):
                return lzma.decompress(payload).decode("utf-8")
            if candidate.endswith(".gz"):
                return gzip.decompress(payload).decode("utf-8")
            return payload.decode("utf-8")
        raise RuntimeError(f"Unable to fetch package index from {base_url}")

    def _container_checks(
        self,
        package_versions_by_dist: dict[str, DebianPackageVersion],
    ) -> list[ContainerAvailabilityCheck]:
        checks: list[ContainerAvailabilityCheck] = []
        for dist, package_version in sorted(package_versions_by_dist.items()):
            normalized = normalize_package_version_tag(str(package_version))
            base_ref = f"{OCI_REGISTRY}/{OCI_REPOSITORY}:{dist}-{normalized}"
            checks.append(
                ContainerAvailabilityCheck(
                    dist=dist,
                    artifact_type="docker",
                    candidate_refs=(base_ref,),
                )
            )
            checks.append(
                ContainerAvailabilityCheck(
                    dist=dist,
                    artifact_type="apptainer",
                    candidate_refs=tuple(f"{base_ref}{suffix}" for suffix in APPTAINER_TAG_SUFFIXES),
                )
            )
        return checks

    def _packaging_target_metadata(self) -> dict[str, dict[str, str]]:
        config_path = self.repo_root / ".github" / "plan-ci.json"
        try:
            payload = json.loads(config_path.read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

        catalog = payload.get("profiles", {}).get("packaging", {}).get("catalog", {})
        metadata_by_dist: dict[str, dict[str, str]] = {}
        for target, row in catalog.items():
            if not isinstance(row, dict):
                continue
            dist = str(row.get("dist") or "").strip()
            if not dist:
                continue
            metadata_by_dist[dist] = {
                "target": str(target),
                "flavor": str(row.get("flavor") or "").strip(),
                "version": str(row.get("version") or "").strip(),
            }
        return metadata_by_dist

    def _ensure_container_available(self, check: ContainerAvailabilityCheck) -> None:
        for candidate_ref in check.candidate_refs:
            command = ["docker", "manifest", "inspect", candidate_ref]
            completed = subprocess.run(
                command,
                cwd=self.repo_root,
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if completed.returncode == 0:
                return
        raise RuntimeError(
            f"{check.artifact_type.capitalize()} artifact is not available: "
            f"{', '.join(check.candidate_refs)}"
        )

    def _package_notes(
        self,
        *,
        package_versions_by_dist: dict[str, DebianPackageVersion],
        package_checks: list[PackageAvailabilityCheck],
        channel: str,
        omitted_dists: list[str],
    ) -> str:
        package_names_by_dist: dict[str, list[str]] = {}
        for check in package_checks:
            package_names_by_dist.setdefault(check.dist, []).append(check.package_name)
        target_metadata = self._packaging_target_metadata()

        lines = [
            "## Packages",
            "",
        ]
        if package_versions_by_dist:
            released_labels = []
            for dist in sorted(package_versions_by_dist):
                metadata = target_metadata.get(dist, {})
                flavor = metadata.get("flavor") or detect_flavor(dist)
                distro_version = metadata.get("version")
                if distro_version:
                    released_labels.append(f"{dist} ({flavor} {distro_version})")
                else:
                    released_labels.append(f"{dist} ({flavor})")
            released = ", ".join(released_labels)
            lines.append(f"- APT packages available for: `{released}`")
            lines.append(f"- Docker images available for: `{released}`")
            lines.append(f"- Apptainer images available for: `{released}`")
        if omitted_dists:
            lines.append(f"- Omitted distros in this release: `{', '.join(omitted_dists)}`")
        for dist, version in sorted(package_versions_by_dist.items()):
            metadata = target_metadata.get(dist, {})
            flavor = metadata.get("flavor") or detect_flavor(dist)
            distro_version = metadata.get("version")
            normalized = normalize_package_version_tag(str(version))
            repo_url = f"http://apt.feelpp.org/{flavor}/{dist}"
            package_names = sorted(set(package_names_by_dist.get(dist, [])))
            target_label = f"{flavor}/{dist}"
            if distro_version:
                target_label = f"{target_label} ({distro_version})"
            lines.extend(
                [
                    "",
                    f"### {target_label}",
                    "",
                    f"- APT package version: `{version}`",
                    f"- APT channel: `{channel}`",
                    f"- APT repository: `{repo_url}`",
                    "",
                    "Install with APT:",
                    "",
                    "```bash",
                    "sudo install -d -m 0755 /etc/apt/keyrings",
                    f"curl -fsSL http://apt.feelpp.org/apt.gpg | sudo gpg --dearmor -o /etc/apt/keyrings/feelpp.gpg",
                    (
                        f"echo 'deb [signed-by=/etc/apt/keyrings/feelpp.gpg] {repo_url} {dist} {channel}' "
                        "| sudo tee /etc/apt/sources.list.d/feelpp.list >/dev/null"
                    ),
                    "sudo apt update",
                    f"sudo apt install {' '.join(package_names)}",
                    "```",
                    "",
                    "Docker image:",
                    "",
                    "```bash",
                    f"docker pull {OCI_REGISTRY}/{OCI_REPOSITORY}:{dist}-{normalized}",
                    "```",
                    "",
                    "Apptainer image:",
                    "",
                    "```bash",
                    f"apptainer pull oras://{OCI_REGISTRY}/{OCI_REPOSITORY}:{dist}-{normalized}-sif",
                    "```",
                ]
            )
        return "\n".join(lines)

    def _publication_notes(
        self,
        *,
        rows: int | None = None,
        since: str | None = None,
    ) -> str:
        service = HalPublicationService(repo_root=self.repo_root)
        try:
            publications = service.fetch(rows=rows, since=since)
        except RuntimeError:
            return ""
        return service.format_markdown(publications)
