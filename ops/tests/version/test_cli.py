from __future__ import annotations

from unittest import mock
import io
import json
import tempfile
import unittest

from feelpp.ops.version.cli import build_parser as build_version_parser, main as version_main

from ops.tests.version.support import VersionRepoMixin


class CliTests(VersionRepoMixin, unittest.TestCase):
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

    def test_cli_parses_release_pretty_flag(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["release", "0.111.0-preview.13", "--pretty"])
        self.assertTrue(args.pretty)

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


if __name__ == "__main__":
    unittest.main()

