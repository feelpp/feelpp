from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import io
import json
import tempfile
import unittest
from unittest import mock

from feelpp.ops.version.cli import build_parser as build_version_parser, main as version_main
from feelpp.ops.version.models import GitHubContributor
from feelpp.ops.version.publications import HalPublication, HalPublicationService
from feelpp.ops.version.release import ReleaseService
from feelpp.ops.version.repository import VersionRepository


class VersionTests(unittest.TestCase):
    def make_repo(self, tmpdir: str) -> Path:
        repo_root = Path(tmpdir) / "repo"
        (repo_root / ".github").mkdir(parents=True)
        (repo_root / "toolboxes" / "cmake").mkdir(parents=True)
        (repo_root / "mor" / "cmake").mkdir(parents=True)
        (repo_root / "packaging" / "manifest").mkdir(parents=True)
        (repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian").mkdir(parents=True)
        (repo_root / "packaging" / "debian" / "feelpp" / "resolute" / "debian").mkdir(parents=True)
        (repo_root / "packaging" / "debian" / "feelpp-toolboxes" / "noble" / "debian").mkdir(parents=True)
        (repo_root / "packaging" / "debian" / "feelpp-mor" / "noble" / "debian").mkdir(parents=True)
        (repo_root / ".git").mkdir()

        version_template = (
            'set(FEELPP_VERSION_MAJOR "{major}")\n'
            'set(FEELPP_VERSION_MINOR "{minor}")\n'
            'set(FEELPP_VERSION_MICRO "{patch}")\n'
            'set(FEELPP_VERSION_PRERELEASE "{prerelease}")\n'
        )
        for path, values in {
            repo_root / "feelpp.version.cmake": (0, 111, 0, "-preview.13"),
            # Legacy per-component version files may still exist in the tree,
            # but fpp-version should ignore them in favor of the repo root.
            repo_root / "toolboxes" / "cmake" / "feelpp.version.cmake": (0, 108, 0, "-beta.1"),
            repo_root / "mor" / "cmake" / "feelpp.version.cmake": (0, 109, 0, "-beta.1"),
        }.items():
            path.write_text(
                version_template.format(
                    major=values[0],
                    minor=values[1],
                    patch=values[2],
                    prerelease=values[3],
                ),
                encoding="utf-8",
            )

        (repo_root / "codemeta.json").write_text(
            json.dumps(
                {
                    "@context": "https://w3id.org/codemeta/3.0",
                    "@type": "SoftwareSourceCode",
                    "name": "Feel++",
                    "version": "v0.111.0-preview.13",
                    "codeRepository": "https://github.com/feelpp/feelpp",
                    "author": [
                        {"@type": "Person", "givenName": "Christophe", "familyName": "Prud'homme"},
                        {"@type": "Person", "givenName": "Vincent", "familyName": "Chabannes"},
                    ],
                    "contributor": [
                        {"@type": "Person", "givenName": "Thomas", "familyName": "Saigre"},
                    ],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        (repo_root / ".zenodo.json").write_text(
            json.dumps(
                {
                    "title": "Feel++",
                    "version": "v0.111.0-preview.13",
                    "upload_type": "software",
                    "license": "lgpl-3.0",
                    "access_right": "open",
                    "creators": [{"name": "Prud'homme, Christophe"}],
                    "contributors": [{"name": "Saigre, Thomas", "affiliation": "IRMA", "type": "Researcher"}],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        (repo_root / "CITATION.cff").write_text(
            "\n".join(
                [
                    "cff-version: 1.2.0",
                    'message: "If you use this software, please cite it as below."',
                    'title: "Feel++"',
                    "version: v0.111.0-preview.13",
                    "type: software",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

        manifest = "\n".join(
            [
                "version = 1",
                'default_components = ["feelpp", "feelpp-toolboxes", "feelpp-mor"]',
                "",
                "[components.feelpp]",
                'distros = ["noble", "resolute"]',
                "dependencies = []",
                'python_packages = ["python3-feelpp"]',
                "publish = true",
                'package_revision = "2"',
                "",
                '[components."feelpp-toolboxes"]',
                'distros = ["noble"]',
                'dependencies = ["feelpp"]',
                'python_packages = ["python3-feelpp-toolboxes"]',
                "publish = true",
                'package_revision = "4"',
                "",
                '[components."feelpp-mor"]',
                'distros = ["noble"]',
                'dependencies = ["feelpp-toolboxes"]',
                'python_packages = ["python3-feelpp-mor"]',
                "publish = true",
                'package_revision = "5"',
                "",
            ]
        )
        (repo_root / "packaging" / "manifest" / "components.toml").write_text(manifest, encoding="utf-8")
        (repo_root / ".github" / "plan-ci.json").write_text(
            json.dumps(
                {
                    "profiles": {
                        "packaging": {
                            "catalog": {
                                "ubuntu:noble": {"flavor": "ubuntu", "dist": "noble", "version": "24.04"},
                                "ubuntu:resolute": {"flavor": "ubuntu", "dist": "resolute", "version": "26.04"},
                                "debian:trixie": {"flavor": "debian", "dist": "trixie", "version": "13"},
                            }
                        }
                    },
                    "releaseNotes": {
                        "publications": {
                            "hal": {
                                "collections": ["FEEL", "CEMOSIS"],
                                "rows": 4,
                            }
                        }
                    },
                }
            ),
            encoding="utf-8",
        )

        changelog_template = (
            "{source} ({version}) unstable; urgency=medium\n\n"
            "  * Existing entry\n\n"
            " -- Test User <test@example.com>  Mon, 01 Jan 2024 00:00:00 +0000\n"
        )
        changelogs = {
            repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog": (
                "feelpp",
                "0.111.0~preview.13-1",
            ),
            repo_root / "packaging" / "debian" / "feelpp" / "resolute" / "debian" / "changelog": (
                "feelpp",
                "0.111.0~preview.13-1",
            ),
            repo_root / "packaging" / "debian" / "feelpp-toolboxes" / "noble" / "debian" / "changelog": (
                "feelpp-toolboxes",
                "0.111.0~preview.13-1",
            ),
            repo_root / "packaging" / "debian" / "feelpp-mor" / "noble" / "debian" / "changelog": (
                "feelpp-mor",
                "0.111.0~preview.13-1",
            ),
        }
        for path, (source, version) in changelogs.items():
            path.write_text(changelog_template.format(source=source, version=version), encoding="utf-8")
        return repo_root

    def test_show_reads_cmake_and_package_versions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            state = VersionRepository(repo_root).read_state()

            self.assertEqual(str(state.canonical_upstream_version()), "0.111.0-preview.13")
            self.assertEqual(state.package_record("feelpp", "noble").version.revision, "2")
            self.assertEqual(state.package_record("feelpp-toolboxes", "noble").version.revision, "4")
            self.assertEqual(len(state.cmake_versions), 1)
            self.assertEqual(len(state.metadata_versions), 3)
            self.assertEqual(state.package_record("feelpp", "noble").origin, "manifest")
            self.assertEqual(state.changelog_record("feelpp", "noble").origin, "changelog")
            self.assertTrue(state.as_dict()["consistency"]["metadata_versions"])

    def test_bump_updates_all_version_targets(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)
            fixed_time = datetime(2026, 4, 8, 12, 0, 0, tzinfo=timezone.utc)

            repository.bump_upstream(
                version=repository.read_state().canonical_upstream_version().parse("0.112.0"),
                message="New upstream release",
                timestamp=fixed_time,
            )

            self.assertIn('set(FEELPP_VERSION_MINOR "112")', (repo_root / "feelpp.version.cmake").read_text(encoding="utf-8"))
            self.assertIn('"version": "v0.112.0"', (repo_root / "codemeta.json").read_text(encoding="utf-8"))
            self.assertIn('"version": "v0.112.0"', (repo_root / ".zenodo.json").read_text(encoding="utf-8"))
            self.assertIn("version: v0.112.0", (repo_root / "CITATION.cff").read_text(encoding="utf-8"))
            manifest_text = (repo_root / "packaging" / "manifest" / "components.toml").read_text(encoding="utf-8")
            self.assertIn('package_revision = "1"', manifest_text)
            self.assertIn("feelpp (0.112.0-1)", (repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0])
            self.assertIn("feelpp-toolboxes (0.112.0-1)", (repo_root / "packaging" / "debian" / "feelpp-toolboxes" / "noble" / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0])

    def test_bump_dry_run_does_not_persist_changes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)
            version_file = repo_root / "feelpp.version.cmake"
            codemeta_file = repo_root / "codemeta.json"
            zenodo_file = repo_root / ".zenodo.json"
            citation_file = repo_root / "CITATION.cff"
            manifest_file = repo_root / "packaging" / "manifest" / "components.toml"
            changelog_file = repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog"
            before = {
                version_file: version_file.read_text(encoding="utf-8"),
                codemeta_file: codemeta_file.read_text(encoding="utf-8"),
                zenodo_file: zenodo_file.read_text(encoding="utf-8"),
                citation_file: citation_file.read_text(encoding="utf-8"),
                manifest_file: manifest_file.read_text(encoding="utf-8"),
                changelog_file: changelog_file.read_text(encoding="utf-8"),
            }

            preview = repository.bump_upstream(
                version=repository.read_state().canonical_upstream_version().parse("0.112.0"),
                dry_run=True,
            )

            self.assertEqual(str(preview.canonical_upstream_version()), "0.112.0")
            self.assertEqual(str(preview.package_record("feelpp", "noble").version), "0.112.0-1")
            for path, expected in before.items():
                self.assertEqual(path.read_text(encoding="utf-8"), expected)

    def test_sync_metadata_updates_release_contributors_only_for_json_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)

            repository.sync_metadata(
                contributors=(
                    GitHubContributor(login="prudhomm", name="Christophe Prud'homme"),
                    GitHubContributor(login="Philand", name="Philippe Pincon"),
                    GitHubContributor(login="t-saigre", name="Thomas Saigre"),
                )
            )

            zenodo_payload = json.loads((repo_root / ".zenodo.json").read_text(encoding="utf-8"))
            self.assertEqual(
                zenodo_payload["contributors"],
                [
                    {"name": "Pincon, Philippe", "type": "Researcher"},
                    {"name": "Saigre, Thomas", "affiliation": "IRMA", "type": "Researcher"},
                ],
            )

            codemeta_payload = json.loads((repo_root / "codemeta.json").read_text(encoding="utf-8"))
            self.assertEqual(
                codemeta_payload["contributor"],
                [
                    {"@type": "Person", "givenName": "Philippe", "familyName": "Pincon"},
                    {"@type": "Person", "givenName": "Thomas", "familyName": "Saigre"},
                ],
            )

            self.assertNotIn("prudhomm", json.dumps(zenodo_payload))
            self.assertNotIn("prudhomm", json.dumps(codemeta_payload))

    def test_revision_bump_increments_debian_revision_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)
            fixed_time = datetime(2026, 4, 8, 12, 0, 0, tzinfo=timezone.utc)

            repository.bump_revision(timestamp=fixed_time)

            manifest_text = (repo_root / "packaging" / "manifest" / "components.toml").read_text(encoding="utf-8")
            self.assertIn('package_revision = "3"', manifest_text)
            self.assertIn('package_revision = "5"', manifest_text)
            self.assertIn('package_revision = "6"', manifest_text)
            self.assertIn("feelpp (0.111.0~preview.13-3)", (repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0])
            self.assertIn('set(FEELPP_VERSION_MINOR "111")', (repo_root / "feelpp.version.cmake").read_text(encoding="utf-8"))

    def test_revision_bump_can_target_specific_dists(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)
            fixed_time = datetime(2026, 4, 8, 12, 0, 0, tzinfo=timezone.utc)

            state = repository.bump_revision(
                timestamp=fixed_time,
                dists=("noble",),
            )

            manifest_text = (repo_root / "packaging" / "manifest" / "components.toml").read_text(encoding="utf-8")
            self.assertIn('package_revision = "2"', manifest_text)
            self.assertIn('package_revision_by_dist = { noble = "3" }', manifest_text)
            self.assertEqual(str(state.package_record("feelpp", "noble").version), "0.111.0~preview.13-3")
            self.assertEqual(str(state.package_record("feelpp", "resolute").version), "0.111.0~preview.13-2")
            self.assertIn(
                "feelpp (0.111.0~preview.13-3)",
                (
                    repo_root
                    / "packaging"
                    / "debian"
                    / "feelpp"
                    / "noble"
                    / "debian"
                    / "changelog"
                ).read_text(encoding="utf-8").splitlines()[0],
            )
            self.assertIn(
                "feelpp (0.111.0~preview.13-1)",
                (
                    repo_root
                    / "packaging"
                    / "debian"
                    / "feelpp"
                    / "resolute"
                    / "debian"
                    / "changelog"
                ).read_text(encoding="utf-8").splitlines()[0],
            )

    def test_revision_bump_dry_run_does_not_persist_changes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)
            manifest_file = repo_root / "packaging" / "manifest" / "components.toml"
            noble_changelog = repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog"
            resolute_changelog = repo_root / "packaging" / "debian" / "feelpp" / "resolute" / "debian" / "changelog"
            before = {
                manifest_file: manifest_file.read_text(encoding="utf-8"),
                noble_changelog: noble_changelog.read_text(encoding="utf-8"),
                resolute_changelog: resolute_changelog.read_text(encoding="utf-8"),
            }

            preview = repository.bump_revision(
                dists=("noble",),
                dry_run=True,
            )

            self.assertEqual(str(preview.package_record("feelpp", "noble").version), "0.111.0~preview.13-3")
            self.assertEqual(str(preview.package_record("feelpp", "resolute").version), "0.111.0~preview.13-2")
            for path, expected in before.items():
                self.assertEqual(path.read_text(encoding="utf-8"), expected)

    def test_sync_changelogs_aligns_headers_with_manifest_versions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)
            fixed_time = datetime(2026, 4, 8, 12, 0, 0, tzinfo=timezone.utc)

            state = repository.read_state()
            self.assertFalse(state.as_dict()["consistency"]["changelog_sync"])

            state = repository.sync_changelogs(timestamp=fixed_time)

            self.assertTrue(state.as_dict()["consistency"]["changelog_sync"])
            self.assertIn(
                "feelpp (0.111.0~preview.13-2)",
                (
                    repo_root
                    / "packaging"
                    / "debian"
                    / "feelpp"
                    / "noble"
                    / "debian"
                    / "changelog"
                ).read_text(encoding="utf-8").splitlines()[0],
            )

    def test_sync_dry_run_does_not_persist_changes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)
            changelog_file = repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog"
            before = changelog_file.read_text(encoding="utf-8")

            state = repository.sync_changelogs(dry_run=True)

            self.assertTrue(state.as_dict()["consistency"]["changelog_sync"])
            self.assertEqual(changelog_file.read_text(encoding="utf-8"), before)

    def test_release_dry_run_builds_plan_without_publishing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            manifest_path = repo_root / "packaging" / "manifest" / "components.toml"
            manifest_path.write_text(
                manifest_path.read_text(encoding="utf-8")
                .replace('package_revision = "4"', 'package_revision = "2"')
                .replace('package_revision = "5"', 'package_revision = "2"'),
                encoding="utf-8",
            )
            service = ReleaseService(repo_root)

            git_outputs = {
                ("branch", "--show-current"): "main\n",
                ("rev-parse", "HEAD"): "abc123\n",
                ("status", "--short"): "",
                ("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"): "origin/main\n",
                ("rev-parse", "origin/main"): "abc123\n",
                ("tag", "--list", "v0.111.0-preview.13"): "",
                ("tag", "--merged", "HEAD", "--list", "v*"): "v0.110.0\nv0.111.0-preview.11\n",
                ("remote", "get-url", "origin"): "https://github.com/feelpp/feelpp.git\n",
            }

            def fake_git_capture(args: list[str], *, check: bool = True) -> str:
                key = tuple(args)
                if key not in git_outputs:
                    raise AssertionError(f"Unexpected git command: {args}")
                return git_outputs[key]

            with mock.patch.object(service, "_git_capture", side_effect=fake_git_capture):
                with mock.patch.object(service, "_ensure_github_checks_green") as gh_checks:
                    with mock.patch.object(service, "_ensure_package_available") as apt_checks:
                        with mock.patch.object(service, "_ensure_container_available") as container_checks:
                            with mock.patch.object(
                                service,
                                "_generated_notes_preview",
                                return_value="## What's Changed\n* Fix packaging",
                            ):
                                with mock.patch.object(
                                    service,
                                    "_publication_notes",
                                    return_value="## Recent Publications using Feel++\n\n- [Paper](https://hal.science/hal-1)",
                                ):
                                    plan = service.execute_release("0.111.0-preview.13", dry_run=True)

            gh_checks.assert_called_once()
            self.assertGreaterEqual(apt_checks.call_count, 1)
            self.assertGreaterEqual(container_checks.call_count, 1)
            self.assertEqual(plan.tag, "v0.111.0-preview.13")
            self.assertIn("## Packages", plan.package_notes)
            self.assertIn(
                "APT packages available for: `noble (ubuntu 24.04), resolute (ubuntu 26.04)`",
                plan.package_notes,
            )
            self.assertIn("### ubuntu/noble (24.04)", plan.package_notes)
            self.assertIn("sudo apt install python3-feelpp", plan.package_notes)
            self.assertIn("## Recent Publications using Feel++", plan.package_notes)
            self.assertIn("docker pull ghcr.io/feelpp/feelpp:noble-v0.111.0-preview.13", plan.package_notes)
            self.assertIn("apptainer pull oras://ghcr.io/feelpp/feelpp:noble-v0.111.0-preview.13-sif", plan.package_notes)
            self.assertIn("## What's Changed", plan.generated_notes_preview)

    def test_release_prepare_rejects_stale_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            (repo_root / "codemeta.json").write_text(
                json.dumps(
                    {
                        "@context": "https://w3id.org/codemeta/3.0",
                        "@type": "SoftwareSourceCode",
                        "name": "Feel++",
                        "version": "v0.111.0-preview.12",
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
            service = ReleaseService(repo_root)

            with mock.patch.object(service, "_git_capture", side_effect=AssertionError("unexpected git")):
                synced = service.repository.sync_metadata()

            self.assertEqual(str(synced.canonical_upstream_version()), "0.111.0-preview.13")
            self.assertIn('"version": "v0.111.0-preview.13"', (repo_root / "codemeta.json").read_text(encoding="utf-8"))

    def test_release_dry_run_can_scope_distros_and_note_omissions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            manifest_path = repo_root / "packaging" / "manifest" / "components.toml"
            manifest_path.write_text(
                manifest_path.read_text(encoding="utf-8")
                .replace('package_revision = "4"', 'package_revision = "2"')
                .replace('package_revision = "5"', 'package_revision = "2"'),
                encoding="utf-8",
            )
            service = ReleaseService(repo_root)

            git_outputs = {
                ("branch", "--show-current"): "main\n",
                ("rev-parse", "HEAD"): "abc123\n",
                ("status", "--short"): "",
                ("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"): "origin/main\n",
                ("rev-parse", "origin/main"): "abc123\n",
                ("tag", "--list", "v0.111.0-preview.13"): "",
                ("tag", "--merged", "HEAD", "--list", "v*"): "v0.110.0\nv0.111.0-preview.11\n",
                ("remote", "get-url", "origin"): "https://github.com/feelpp/feelpp.git\n",
            }

            def fake_git_capture(args: list[str], *, check: bool = True) -> str:
                key = tuple(args)
                if key not in git_outputs:
                    raise AssertionError(f"Unexpected git command: {args}")
                return git_outputs[key]

            with mock.patch.object(service, "_git_capture", side_effect=fake_git_capture):
                with mock.patch.object(service, "_ensure_github_checks_green") as gh_checks:
                    with mock.patch.object(service, "_ensure_package_available") as apt_checks:
                        with mock.patch.object(service, "_ensure_container_available") as container_checks:
                            with mock.patch.object(
                                service,
                                "_generated_notes_preview",
                                return_value="## What's Changed\n* Scoped release",
                            ):
                                with mock.patch.object(service, "_publication_notes", return_value=""):
                                    plan = service.execute_release(
                                        "0.111.0-preview.13",
                                        dry_run=True,
                                        dists=("noble",),
                                    )

            gh_checks.assert_called_once()
            self.assertEqual(apt_checks.call_count, 3)
            self.assertEqual(container_checks.call_count, 2)
            self.assertEqual({check.dist for check in plan.package_checks}, {"noble"})
            self.assertEqual({check.dist for check in plan.container_checks}, {"noble"})
            self.assertIn("APT packages available for: `noble (ubuntu 24.04)`", plan.package_notes)
            self.assertIn("Omitted distros in this release: `resolute`", plan.package_notes)

    def test_generated_notes_preview_uses_github_release_notes_api(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = ReleaseService(repo_root)
            with mock.patch(
                "feelpp.ops.version.release.run_capture",
                return_value=json.dumps({"body": "## What's Changed\n* closes #1"}),
            ) as run_capture_mock:
                preview = service._generated_notes_preview(
                    repo_slug="feelpp/feelpp",
                    tag="v0.111.0-preview.13",
                    head_sha="abc123",
                    previous_tag="v0.111.0-preview.12",
                )

        self.assertIn("## What's Changed", preview)
        self.assertEqual(
            run_capture_mock.call_args.args[0],
            [
                "gh",
                "api",
                "repos/feelpp/feelpp/releases/generate-notes",
                "-X",
                "POST",
                "-f",
                "tag_name=v0.111.0-preview.13",
                "-f",
                "target_commitish=abc123",
                "-f",
                "previous_tag_name=v0.111.0-preview.12",
            ],
        )

    def test_previous_tag_for_prerelease_uses_previous_prerelease_in_same_series(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = ReleaseService(repo_root)

            with mock.patch.object(
                service,
                "_git_capture",
                return_value="v0.110.0\nv0.111.0-preview.1\nv0.111.0-preview.11\n",
            ):
                previous = service._previous_tag("v0.111.0-preview.13")

        self.assertEqual(previous, "v0.111.0-preview.11")

    def test_previous_tag_for_stable_release_skips_prereleases(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = ReleaseService(repo_root)

            with mock.patch.object(
                service,
                "_git_capture",
                return_value="v0.110.0\nv0.111.0-preview.11\n",
            ):
                previous = service._previous_tag("v0.111.0")

        self.assertEqual(previous, "v0.110.0")

    def test_previous_tag_for_patch_release_uses_previous_stable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = ReleaseService(repo_root)

            with mock.patch.object(
                service,
                "_git_capture",
                return_value="v0.110.0\nv0.111.0-preview.11\nv0.111.0\n",
            ):
                previous = service._previous_tag("v0.111.1")

        self.assertEqual(previous, "v0.111.0")

    def test_release_contributors_are_derived_from_compare_range(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = ReleaseService(repo_root)
            plan = mock.Mock(
                previous_tag="v0.111.0-preview.12",
                repo_slug="feelpp/feelpp",
                head_sha="abc123",
            )

            responses = [
                {
                    "commits": [
                        {"author": {"login": "prudhomm"}},
                        {"author": {"login": "Philand"}},
                        {"author": {"login": "dependabot[bot]"}},
                        {"author": {"login": "Philand"}},
                        {"author": None},
                    ]
                },
                {"login": "prudhomm", "name": "Christophe Prud'homme", "type": "User"},
                {"login": "Philand", "name": "Philippe Pincon", "type": "User"},
                {"login": "dependabot[bot]", "name": "dependabot", "type": "Bot"},
            ]

            def fake_run_capture(args, *, cwd=None, check=True):
                return json.dumps(responses.pop(0))

            with mock.patch("feelpp.ops.version.release.run_capture", side_effect=fake_run_capture):
                contributors = service._release_contributors(plan)

            self.assertEqual(
                contributors,
                (
                    GitHubContributor(login="prudhomm", name="Christophe Prud'homme"),
                    GitHubContributor(login="Philand", name="Philippe Pincon"),
                ),
            )

    def test_release_contributors_require_previous_tag(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = ReleaseService(repo_root)
            plan = mock.Mock(previous_tag=None, repo_slug="feelpp/feelpp", head_sha="abc123")

            contributors = service._release_contributors(plan)

            self.assertIsNone(contributors)

    def test_execute_release_commits_metadata_before_tag(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            manifest_path = repo_root / "packaging" / "manifest" / "components.toml"
            manifest_path.write_text(
                manifest_path.read_text(encoding="utf-8")
                .replace('package_revision = "4"', 'package_revision = "2"')
                .replace('package_revision = "5"', 'package_revision = "2"'),
                encoding="utf-8",
            )
            (repo_root / "codemeta.json").write_text(
                (repo_root / "codemeta.json").read_text(encoding="utf-8").replace("v0.111.0-preview.13", "v0.111.0-preview.12"),
                encoding="utf-8",
            )
            service = ReleaseService(repo_root)

            git_outputs = {
                ("branch", "--show-current"): "main\n",
                ("rev-parse", "HEAD"): "def456\n",
                ("status", "--short"): "",
                ("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"): "origin/main\n",
                ("rev-parse", "origin/main"): "def456\n",
                ("tag", "--list", "v0.111.0-preview.13"): "",
                ("tag", "--merged", "HEAD", "--list", "v*"): "v0.110.0\nv0.111.0-preview.11\n",
                ("remote", "get-url", "origin"): "https://github.com/feelpp/feelpp.git\n",
                ("status", "--short", "--", "codemeta.json", ".zenodo.json", "CITATION.cff"): " M codemeta.json\n",
            }
            head_calls = {"count": 0}

            def fake_git_capture(args: list[str], *, check: bool = True) -> str:
                key = tuple(args)
                if key == ("rev-parse", "HEAD"):
                    head_calls["count"] += 1
                    return "def456\n" if head_calls["count"] <= 2 else "fedcba\n"
                if key not in git_outputs:
                    raise AssertionError(f"Unexpected git command: {args}")
                return git_outputs[key]

            checked_commands: list[list[str]] = []

            def fake_run_checked(args, *, cwd=None, **kwargs) -> None:
                checked_commands.append([str(arg) for arg in args])

            with mock.patch.object(service, "_git_capture", side_effect=fake_git_capture):
                with mock.patch.object(service, "_ensure_github_checks_green"):
                    with mock.patch.object(service, "_ensure_package_available"):
                        with mock.patch.object(service, "_ensure_container_available"):
                            with mock.patch.object(service, "_generated_notes_preview", return_value="## What's Changed\n* Fix packaging"):
                                with mock.patch.object(service, "_publication_notes", return_value=""):
                                    with mock.patch.object(
                                        service,
                                        "_release_contributors",
                                        return_value=(GitHubContributor(login="Philand", name="Philippe Pincon"),),
                                    ):
                                        with mock.patch("feelpp.ops.version.release.run_checked", side_effect=fake_run_checked):
                                            plan = service.execute_release("0.111.0-preview.13", dry_run=False)

            self.assertEqual(plan.head_sha, "fedcba")
            self.assertIn(["git", "add", "--", "codemeta.json", ".zenodo.json", "CITATION.cff"], checked_commands)
            self.assertIn(
                ["git", "commit", "-m", "chore(release): sync metadata for v0.111.0-preview.13 [ci skip]"],
                checked_commands,
            )
            self.assertIn(["git", "push", "origin", "main"], checked_commands)
            self.assertIn(["git", "tag", "-a", "v0.111.0-preview.13", "-m", "v0.111.0-preview.13", "fedcba"], checked_commands)
            self.assertIn('"version": "v0.111.0-preview.13"', (repo_root / "codemeta.json").read_text(encoding="utf-8"))
            self.assertIn('"familyName": "Pincon"', (repo_root / "codemeta.json").read_text(encoding="utf-8"))
            self.assertIn('"name": "Pincon, Philippe"', (repo_root / ".zenodo.json").read_text(encoding="utf-8"))

    def test_hal_publication_service_uses_feel_and_cemosis_collections(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = HalPublicationService(repo_root)
            payload = {
                "response": {
                    "docs": [
                        {
                            "docid": "1",
                            "title_s": ["A Feel++ paper"],
                            "uri_s": "https://hal.science/hal-1",
                            "producedDate_tdate": "2026-01-28T00:00:00Z",
                            "producedDateY_i": 2026,
                            "docType_s": "UNDEFINED",
                            "authFullName_s": ["Alice", "Bob", "Carol", "Dave"],
                            "journalTitle_s": "HAL Preprint",
                        }
                    ]
                }
            }

            class FakeResponse:
                def __enter__(self):
                    return self

                def __exit__(self, exc_type, exc, tb):
                    return False

                def read(self):
                    return json.dumps(payload).encode("utf-8")

            with mock.patch("feelpp.ops.version.publications.urlopen", return_value=FakeResponse()) as urlopen_mock:
                publications = service.fetch()

            request = urlopen_mock.call_args.args[0]
            self.assertIn("collCode_s%3A%28FEEL+OR+CEMOSIS%29", request.full_url)
            self.assertEqual(len(publications), 1)
            self.assertEqual(publications[0].title, "A Feel++ paper")
            markdown = service.format_markdown(publications)
            self.assertIn("## Recent Publications using Feel++", markdown)
            self.assertIn("Alice, Bob, Carol, et al.", markdown)
            self.assertIn("Preprint", markdown)

    def test_cli_exposes_revision_bump_command(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["revision", "bump"])
        self.assertEqual(args.revision_command, "bump")
        self.assertEqual(args.dist, [])
        self.assertFalse(args.dry_run)

    def test_cli_exposes_sync_command(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["sync"])
        self.assertEqual(args.command, "sync")
        self.assertEqual(args.dist, [])
        self.assertFalse(args.dry_run)

    def test_cli_parses_revision_bump_dist_and_dry_run(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["revision", "bump", "--dist", "noble", "--dist", "trixie", "--dry-run"])
        self.assertEqual(args.dist, ["noble", "trixie"])
        self.assertTrue(args.dry_run)

    def test_cli_parses_bump_dry_run(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["bump", "0.112.0", "--dry-run"])
        self.assertTrue(args.dry_run)

    def test_cli_parses_release_dist_and_dry_run(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["release", "0.111.0-preview.13", "--dist", "noble", "--dist", "trixie", "--dry-run"])
        self.assertEqual(args.dist, ["noble", "trixie"])
        self.assertTrue(args.dry_run)
        self.assertFalse(args.pretty)

    def test_cli_parses_release_publication_flags(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(
            [
                "release",
                "0.111.0-preview.13",
                "--publications-rows",
                "7",
                "--publications-since",
                "2026-01-01",
            ]
        )
        self.assertEqual(args.publications_rows, 7)
        self.assertEqual(args.publications_since, "2026-01-01")

    def test_cli_parses_release_pretty_flag(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["release", "0.111.0-preview.13", "--pretty"])
        self.assertTrue(args.pretty)

    def test_cli_parses_publications_pretty_flag(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["publications", "--pretty", "--collection", "FEEL"])
        self.assertTrue(args.pretty)
        self.assertEqual(args.collection, ["FEEL"])

    def test_main_prints_error_for_invalid_semver(self) -> None:
        stderr = io.StringIO()
        with mock.patch("sys.stderr", stderr):
            rc = version_main(["bump", "invalid-semver"])

        self.assertEqual(rc, 1)
        self.assertIn("Invalid semantic version", stderr.getvalue())

    def test_show_command_prints_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            stdout = io.StringIO()
            with mock.patch("sys.stdout", stdout):
                rc = version_main(["--repo-root", str(repo_root), "show"])

            self.assertEqual(rc, 0)
            payload = json.loads(stdout.getvalue())
            self.assertEqual(payload["consistency"]["cmake_versions"], True)
            self.assertEqual(payload["package_versions"][0]["origin"], "manifest")

    def test_release_command_pretty_prints_plan(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            stdout = io.StringIO()

            class FakePlan:
                title = "v0.111.0-preview.13"
                tag = "v0.111.0-preview.13"
                release_kind = "upstream"
                branch = "develop"
                channel = "latest"
                head_sha = "abc123"
                previous_tag = "v0.111.0-preview.12"
                prerelease = True
                dry_run = True
                package_notes = "## Packages\n\n- APT packages available for: `noble, trixie`"
                generated_notes_preview = "* abc123 Test commit"
                package_checks = (
                    mock.Mock(
                        dist="noble",
                        package_name="python3-feelpp",
                        expected_version="0.111.0~preview.13-2",
                        url="http://apt.feelpp.org/ubuntu/noble",
                    ),
                )
                container_checks = (
                    mock.Mock(
                        dist="noble",
                        artifact_type="docker",
                        candidate_refs=("ghcr.io/feelpp/feelpp:noble-v0.111.0-preview.13",),
                    ),
                )

            with mock.patch("feelpp.ops.version.cli.ReleaseService") as service_cls:
                service_cls.return_value.execute_release.return_value = FakePlan()
                with mock.patch("sys.stdout", stdout):
                    rc = version_main(["--repo-root", str(repo_root), "release", "0.111.0-preview.13", "--pretty"])

            self.assertEqual(rc, 0)
            output = stdout.getvalue()
            self.assertIn("Release: v0.111.0-preview.13", output)
            self.assertIn("Package notes:", output)
            self.assertIn("Generated notes preview:", output)
            self.assertNotIn("Package checks:", output)
            self.assertNotIn("Container checks:", output)

    def test_publications_command_pretty_prints_markdown(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            stdout = io.StringIO()
            publication = HalPublication(
                docid="1",
                title="A Feel++ paper",
                url="https://hal.science/hal-1",
                produced_at="2026-01-28T00:00:00Z",
                year=2026,
                doc_type="UNDEFINED",
                venue="HAL Preprint",
                authors=("Alice", "Bob"),
                doi=None,
            )

            with mock.patch("feelpp.ops.version.cli.HalPublicationService") as service_cls:
                service_cls.return_value.fetch.return_value = [publication]
                service_cls.return_value.format_markdown.return_value = "## Recent Publications using Feel++\n\n- [A Feel++ paper](https://hal.science/hal-1)"
                with mock.patch("sys.stdout", stdout):
                    rc = version_main(["--repo-root", str(repo_root), "publications", "--pretty"])

            self.assertEqual(rc, 0)
            self.assertIn("## Recent Publications using Feel++", stdout.getvalue())

    def test_release_command_passes_publication_overrides(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            stdout = io.StringIO()

            class FakePlan:
                title = "v0.111.0-preview.13"
                tag = "v0.111.0-preview.13"
                release_kind = "upstream"
                branch = "develop"
                channel = "latest"
                head_sha = "abc123"
                previous_tag = "v0.111.0-preview.12"
                prerelease = True
                dry_run = True
                package_notes = "## Packages"
                generated_notes_preview = "## What's Changed"

            with mock.patch("feelpp.ops.version.cli.ReleaseService") as service_cls:
                service_cls.return_value.execute_release.return_value = FakePlan()
                with mock.patch("sys.stdout", stdout):
                    rc = version_main(
                        [
                            "--repo-root",
                            str(repo_root),
                            "release",
                            "0.111.0-preview.13",
                            "--publications-rows",
                            "7",
                            "--publications-since",
                            "2026-01-01",
                            "--pretty",
                        ]
                    )

            self.assertEqual(rc, 0)
            service_cls.return_value.execute_release.assert_called_once_with(
                "0.111.0-preview.13",
                dry_run=False,
                dists=None,
                publication_rows=7,
                publication_since="2026-01-01",
            )


if __name__ == "__main__":
    unittest.main()
