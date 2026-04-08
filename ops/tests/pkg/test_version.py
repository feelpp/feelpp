from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import io
import json
import tempfile
import unittest
from unittest import mock

from feelpp.ops.version.cli import build_parser as build_version_parser, main as version_main
from feelpp.ops.version.release import ReleaseService
from feelpp.ops.version.repository import VersionRepository


class VersionTests(unittest.TestCase):
    def make_repo(self, tmpdir: str) -> Path:
        repo_root = Path(tmpdir) / "repo"
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
            repo_root / "toolboxes" / "cmake" / "feelpp.version.cmake": (0, 111, 0, "-preview.13"),
            repo_root / "mor" / "cmake" / "feelpp.version.cmake": (0, 111, 0, "-preview.13"),
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
                "",
                '[components."feelpp-toolboxes"]',
                'distros = ["noble"]',
                'dependencies = ["feelpp"]',
                'python_packages = ["python3-feelpp-toolboxes"]',
                "publish = true",
                "",
                '[components."feelpp-mor"]',
                'distros = ["noble"]',
                'dependencies = ["feelpp-toolboxes"]',
                'python_packages = ["python3-feelpp-mor"]',
                "publish = true",
                "",
            ]
        )
        (repo_root / "packaging" / "manifest" / "components.toml").write_text(manifest, encoding="utf-8")

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
            self.assertEqual(state.package_record("feelpp", "noble").version.revision, "1")

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
            self.assertIn("feelpp (0.112.0-1)", (repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0])
            self.assertIn("feelpp-toolboxes (0.112.0-1)", (repo_root / "packaging" / "debian" / "feelpp-toolboxes" / "noble" / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0])

    def test_revision_bump_increments_debian_revision_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            repository = VersionRepository(repo_root)
            fixed_time = datetime(2026, 4, 8, 12, 0, 0, tzinfo=timezone.utc)

            repository.bump_revision(timestamp=fixed_time)

            self.assertIn("feelpp (0.111.0~preview.13-2)", (repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog").read_text(encoding="utf-8").splitlines()[0])
            self.assertIn('set(FEELPP_VERSION_MINOR "111")', (repo_root / "feelpp.version.cmake").read_text(encoding="utf-8"))

    def test_release_dry_run_builds_plan_without_publishing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = ReleaseService(repo_root)

            git_outputs = {
                ("branch", "--show-current"): "main\n",
                ("rev-parse", "HEAD"): "abc123\n",
                ("status", "--short"): "",
                ("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"): "origin/main\n",
                ("rev-parse", "origin/main"): "abc123\n",
                ("tag", "--list", "v0.111.0-preview.13"): "",
                ("describe", "--tags", "--abbrev=0", "--match", "v*"): "v0.111.0-preview.12\n",
                ("remote", "get-url", "origin"): "https://github.com/feelpp/feelpp.git\n",
                ("log", "--pretty=format:* %h %s", "v0.111.0-preview.12..HEAD"): "* abc123 Test commit\n",
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
                            plan = service.execute_release("0.111.0-preview.13", dry_run=True)

            gh_checks.assert_called_once()
            self.assertGreaterEqual(apt_checks.call_count, 1)
            self.assertGreaterEqual(container_checks.call_count, 1)
            self.assertEqual(plan.tag, "v0.111.0-preview.13")
            self.assertIn("## Packages", plan.package_notes)
            self.assertIn("Test commit", plan.generated_notes_preview)

    def test_cli_exposes_revision_bump_command(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["revision", "bump"])
        self.assertEqual(args.revision_command, "bump")

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


if __name__ == "__main__":
    unittest.main()
