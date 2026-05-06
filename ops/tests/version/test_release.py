from __future__ import annotations

from unittest import mock
import json
import unittest
import tempfile

from feelpp.ops.version.models import GitHubContributor
from feelpp.ops.version.release import GENERATED_NOTES_UNAVAILABLE, ReleaseService

from ops.tests.version.support import VersionRepoMixin


class ReleaseTests(VersionRepoMixin, unittest.TestCase):
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

    def test_release_contributors_are_derived_from_generated_notes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = ReleaseService(repo_root)
            plan = mock.Mock(
                generated_notes_preview=(
                    "## What's Changed\n"
                    "* Fix packaging by @prudhomm in https://github.com/feelpp/feelpp/pull/1\n"
                    "* Improve release notes by @Philand in https://github.com/feelpp/feelpp/pull/2\n"
                    "## New Contributors\n"
                    "* @Philand made their first contribution in https://github.com/feelpp/feelpp/pull/2\n"
                    "* @dependabot[bot] updated dependencies in https://github.com/feelpp/feelpp/pull/3\n"
                )
            )

            responses = [
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

    def test_release_contributors_are_empty_when_notes_are_unavailable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = ReleaseService(repo_root)
            plan = mock.Mock(generated_notes_preview=GENERATED_NOTES_UNAVAILABLE)

            contributors = service._release_contributors(plan)

            self.assertEqual(contributors, ())

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


if __name__ == "__main__":
    unittest.main()
