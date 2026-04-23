from __future__ import annotations

import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.cli import build_parser, main


class CliTests(unittest.TestCase):
    def repo_root(self) -> Path:
        return Path(__file__).resolve().parents[3]

    def test_inspect_plan_defaults_to_host_engine(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["inspect", "plan", "--dist", "trixie"])
        self.assertEqual(args.engine, "host")

    def test_job_init_defaults_to_host_engine(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["job", "init", "--dist", "noble"])
        self.assertEqual(args.engine, "host")

    def test_publish_cleanup_command_is_available(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["publish", "cleanup", "--dist", "noble"])
        self.assertEqual(args.publish_command, "cleanup")

    def test_spack_env_list_command_is_available(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["spack", "env", "list"])
        self.assertEqual(args.command, "spack")
        self.assertEqual(args.spack_command, "env")
        self.assertEqual(args.spack_env_command, "list")

    def test_spack_image_targets_command_is_available(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["spack", "image", "targets"])
        self.assertEqual(args.command, "spack")
        self.assertEqual(args.spack_command, "image")
        self.assertEqual(args.spack_image_command, "targets")

    def test_top_level_image_targets_command_is_available(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["image", "targets"])
        self.assertEqual(args.command, "image")
        self.assertEqual(args.image_command, "targets")

    def test_top_level_image_bake_command_is_available(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["image", "bake", "--target", "ubuntu:noble"])
        self.assertEqual(args.command, "image")
        self.assertEqual(args.image_command, "bake")
        self.assertEqual(args.target, "ubuntu:noble")

    def test_debian_backend_job_init_command_is_available(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["debian", "job", "init", "--dist", "noble"])
        self.assertEqual(args.command, "debian")
        self.assertEqual(args.job_command, "init")
        self.assertEqual(args.engine, "host")

    def test_fpp_spack_alias_prepends_spack_backend(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest = Path(tmpdir) / "packaging" / "spack" / "environments" / "demo" / "spack.yaml"
            manifest.parent.mkdir(parents=True)
            manifest.write_text("spack:\n  specs: []\n", encoding="utf-8")

            stdout = io.StringIO()
            with mock.patch("sys.argv", ["fpp-spack", "env", "list", "--repo-root", tmpdir]):
                with mock.patch("sys.stdout", stdout):
                    rc = main()

        self.assertEqual(rc, 0)
        self.assertIn('"name": "demo"', stdout.getvalue())

    def test_spack_env_list_preserves_nested_environment_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest = (
                Path(tmpdir)
                / "packaging"
                / "spack"
                / "environments"
                / "cpu"
                / "openmpi"
                / "spack.yaml"
            )
            manifest.parent.mkdir(parents=True)
            manifest.write_text("spack:\n  specs: []\n", encoding="utf-8")

            stdout = io.StringIO()
            with mock.patch("sys.stdout", stdout):
                rc = main(["spack", "env", "list", "--repo-root", tmpdir])

        self.assertEqual(rc, 0)
        self.assertIn('"name": "cpu/openmpi"', stdout.getvalue())
        self.assertIn('"status": "supported"', stdout.getvalue())

    def test_spack_image_targets_reflect_plan_ci_catalog(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            plan_path = Path(tmpdir) / ".github" / "plan-ci.json"
            plan_path.parent.mkdir(parents=True)
            plan_path.write_text(
                """
{
  "profiles": {
    "images": {
      "catalog": {
        "spack:openmpi": {
          "flavor": "spack",
          "dist": "openmpi",
          "version": "latest",
          "docker": "true",
          "image_backend": "spack",
          "image_strategy": "full",
          "base_image": "ubuntu:24.04",
          "oci_dist": "spack-openmpi",
          "spack_environment": "cpu/openmpi"
        },
        "ubuntu:noble": {
          "flavor": "ubuntu",
          "dist": "noble",
          "version": "24.04",
          "docker": "true",
          "image_backend": "apt",
          "image_strategy": "components",
          "base_image": "ubuntu:24.04",
          "oci_dist": "ubuntu-24.04"
        }
      }
    }
  }
}
""".strip()
                + "\n",
                encoding="utf-8",
            )
            manifest = (
                Path(tmpdir)
                / "packaging"
                / "spack"
                / "environments"
                / "cpu"
                / "openmpi"
                / "spack.yaml"
            )
            manifest.parent.mkdir(parents=True)
            manifest.write_text("spack:\n  specs: []\n", encoding="utf-8")

            stdout = io.StringIO()
            with mock.patch("sys.stdout", stdout):
                rc = main(["spack", "image", "targets", "--repo-root", tmpdir])

        self.assertEqual(rc, 0)
        self.assertIn('"target": "spack:openmpi"', stdout.getvalue())
        self.assertIn('"environment": "cpu/openmpi"', stdout.getvalue())
        self.assertIn('"supported": true', stdout.getvalue())

    def test_spack_image_bake_writes_bake_ready_context(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            plan_path = repo_root / ".github" / "plan-ci.json"
            plan_path.parent.mkdir(parents=True)
            plan_path.write_text(
                """
{
  "profiles": {
    "images": {
      "catalog": {
        "spack:openmpi": {
          "flavor": "spack",
          "dist": "openmpi",
          "version": "latest",
          "docker": "true",
          "image_backend": "spack",
          "image_strategy": "full",
          "base_image": "ubuntu:24.04",
          "oci_dist": "spack-openmpi",
          "spack_environment": "cpu/openmpi"
        }
      }
    }
  }
}
""".strip()
                + "\n",
                encoding="utf-8",
            )
            manifest = repo_root / "packaging" / "spack" / "environments" / "cpu" / "openmpi" / "spack.yaml"
            manifest.parent.mkdir(parents=True)
            manifest.write_text("spack:\n  specs: []\n", encoding="utf-8")
            (repo_root / "packaging" / "spack" / "README.md").write_text("spack docs\n", encoding="utf-8")
            generated_state = manifest.parent / ".spack-env"
            generated_state.mkdir()
            (generated_state / "view-marker.txt").write_text("generated\n", encoding="utf-8")

            stdout = io.StringIO()
            with mock.patch("sys.stdout", stdout):
                rc = main(
                    [
                        "spack",
                        "image",
                        "bake",
                        "--repo-root",
                        tmpdir,
                        "--job-root",
                        str(repo_root / "job"),
                        "--target",
                        "spack:openmpi",
                    ]
                )

            context_dir = repo_root / "job" / "images" / "spack-openmpi"
            dockerfile = context_dir / "Dockerfile"
            bake_file = context_dir / "docker-bake.json"
            self.assertEqual(rc, 0)
            self.assertTrue(dockerfile.is_file())
            self.assertTrue(bake_file.is_file())
            self.assertIn(
                "COPY packaging/spack /opt/feelpp/packaging/spack",
                dockerfile.read_text(encoding="utf-8"),
            )
            self.assertIn('"spack-openmpi"', bake_file.read_text(encoding="utf-8"))
            self.assertIn('"ghcr.io/feelpp/feelpp:spack-openmpi-full-dev"', bake_file.read_text(encoding="utf-8"))
            self.assertIn('"packaging_target": "spack:openmpi"', stdout.getvalue())
            self.assertFalse((context_dir / "packaging" / "spack" / "environments" / "cpu" / "openmpi" / ".spack-env").exists())

    def test_top_level_image_targets_lists_repo_owned_images_profile(self) -> None:
        stdout = io.StringIO()
        with mock.patch("sys.stdout", stdout):
            rc = main(["image", "targets", "--repo-root", str(self.repo_root()), "--backend", "apt"])

        self.assertEqual(rc, 0)
        self.assertIn('"target": "ubuntu:noble"', stdout.getvalue())
        self.assertIn('"image_backend": "apt"', stdout.getvalue())

    def test_top_level_image_bake_writes_component_bake_graph(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.repo_root()
            stdout = io.StringIO()
            with mock.patch("sys.stdout", stdout):
                rc = main(
                    [
                        "image",
                        "bake",
                        "--repo-root",
                        str(repo_root),
                        "--job-root",
                        tmpdir,
                        "--target",
                        "ubuntu:noble",
                    ]
                )

            context_dir = Path(tmpdir) / "images" / "ubuntu-noble"
            env_dockerfile = context_dir / "feelpp-env" / "Dockerfile"
            bake_file = context_dir / "docker-bake.json"
            self.assertEqual(rc, 0)
            self.assertTrue(env_dockerfile.is_file())
            self.assertTrue(bake_file.is_file())
            bake_payload = bake_file.read_text(encoding="utf-8")
            self.assertIn('"feelpp-env"', bake_payload)
            self.assertIn('"toolboxes-runtime"', bake_payload)
            self.assertIn('"full-all"', bake_payload)
            self.assertIn('"feelpp_env_image": "target:feelpp-env"', bake_payload)
            self.assertIn('"docker_bake_command": "docker buildx bake -f ', stdout.getvalue())
            self.assertIn('"available_groups": [', stdout.getvalue())

    def test_main_prints_clean_error_for_expected_packaging_failures(self) -> None:
        stderr = io.StringIO()
        with mock.patch("feelpp.pkg.backends.debian.commands.inspect.context_from_args", return_value=object()):
            with mock.patch("feelpp.pkg.backends.debian.commands.inspect.load_plan", side_effect=ValueError("boom")):
                with mock.patch("sys.stderr", stderr):
                    rc = main(["inspect", "plan", "--dist", "trixie"])

        self.assertEqual(rc, 1)
        self.assertEqual(stderr.getvalue().strip(), "error: boom")


if __name__ == "__main__":
    unittest.main()
