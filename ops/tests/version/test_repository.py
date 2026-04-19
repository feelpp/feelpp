from __future__ import annotations

from datetime import datetime, timezone
import json
import tempfile
import unittest

from feelpp.ops.version.models import GitHubContributor
from feelpp.ops.version.repository import VersionRepository

from ops.tests.version.support import VersionRepoMixin


class RepositoryTests(VersionRepoMixin, unittest.TestCase):
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
            self.assertIn(
                "feelpp (0.112.0-1)",
                (repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0],
            )
            self.assertIn(
                "feelpp-toolboxes (0.112.0-1)",
                (repo_root / "packaging" / "debian" / "feelpp-toolboxes" / "noble" / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0],
            )

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

    def test_sync_metadata_can_clear_release_contributors(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)

            repository.sync_metadata(contributors=())

            zenodo_payload = json.loads((repo_root / ".zenodo.json").read_text(encoding="utf-8"))
            self.assertEqual(zenodo_payload["contributors"], [])

            codemeta_payload = json.loads((repo_root / "codemeta.json").read_text(encoding="utf-8"))
            self.assertEqual(codemeta_payload["contributor"], [])

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
            self.assertIn(
                "feelpp (0.111.0~preview.13-3)",
                (repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0],
            )
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


if __name__ == "__main__":
    unittest.main()

