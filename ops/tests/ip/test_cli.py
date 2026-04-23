from __future__ import annotations

from pathlib import Path
from unittest import mock
import io
import json
import tempfile
import unittest

from feelpp.ops.ip.cli import build_parser, main as ip_main
from feelpp.ops.ip.exporter import build_public_bundle, render_public_bundle
from feelpp.ops.ip.repository import PublicMetadataRepository
from feelpp.ops.ip.yamlio import load_yaml

from ops.tests.version.support import VersionRepoMixin

from .support import write_public_metadata


class IpCliTests(VersionRepoMixin, unittest.TestCase):
    def test_cli_exposes_public_commands(self) -> None:
        parser = build_parser()
        self.assertEqual(parser.parse_args(["validate-public"]).command, "validate-public")
        self.assertEqual(parser.parse_args(["show-public"]).format, "yaml")
        export_args = parser.parse_args(["export-public", "--format", "json", "--out", "bundle.json"])
        self.assertEqual(export_args.format, "json")
        self.assertEqual(export_args.out, "bundle.json")
        self.assertFalse(parser.parse_args(["stats"]).write_metadata)

    def test_validate_public_command_prints_ok(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            write_public_metadata(repo_root)
            stdout = io.StringIO()

            with mock.patch("sys.stdout", stdout):
                rc = ip_main(["--repo-root", str(repo_root), "validate-public"])

            self.assertEqual(rc, 0)
            self.assertEqual(stdout.getvalue().strip(), "metadata/software.public.yml: ok")

    def test_show_public_json_outputs_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            write_public_metadata(repo_root)
            stdout = io.StringIO()

            with mock.patch("sys.stdout", stdout):
                rc = ip_main(["--repo-root", str(repo_root), "show-public", "--format", "json"])

            self.assertEqual(rc, 0)
            payload = json.loads(stdout.getvalue())
            self.assertEqual(payload["software"]["name"], "Feel++")

    def test_export_public_json_writes_release_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            write_public_metadata(repo_root)
            out_path = Path(tmpdir) / "bundle.json"

            rc = ip_main(
                [
                    "--repo-root",
                    str(repo_root),
                    "export-public",
                    "--format",
                    "json",
                    "--out",
                    str(out_path),
                ]
            )

            self.assertEqual(rc, 0)
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["kind"], "feelpp-public-app-bundle")
            self.assertEqual(payload["release"]["version"], "0.111.0-preview.13")
            self.assertEqual(payload["source"]["metadata_path"], "metadata/software.public.yml")

    def test_public_export_is_stable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            write_public_metadata(repo_root)
            repository = PublicMetadataRepository(repo_root=repo_root)

            first = render_public_bundle(build_public_bundle(repository), "yaml")
            second = render_public_bundle(build_public_bundle(repository), "yaml")

            self.assertEqual(first, second)
            self.assertNotIn("generated_at", first)

    def test_stats_write_metadata_updates_metrics_block(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            metadata_path = write_public_metadata(repo_root)

            rc = ip_main(["--repo-root", str(repo_root), "stats", "--write-metadata"])

            self.assertEqual(rc, 0)
            payload = load_yaml(metadata_path)
            self.assertGreater(payload["metrics"]["counted_files"], 0)
            self.assertGreater(payload["metrics"]["approximate_bytes"], 0)


if __name__ == "__main__":
    unittest.main()
