from __future__ import annotations

import io
import json
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

    def test_spack_image_build_accepts_parallelism_overrides(self) -> None:
        parser = build_parser()
        args = parser.parse_args(
            [
                "spack",
                "image",
                "build",
                "--target",
                "spack:openmpi",
                "--spack-build-jobs",
                "24",
                "--spack-concurrent-packages",
                "3",
                "--dry-run",
            ]
        )
        self.assertEqual(args.spack_build_jobs, 24)
        self.assertEqual(args.spack_concurrent_packages, 3)

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

    def test_top_level_image_build_command_is_available(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["image", "build", "--target", "ubuntu:noble", "--push"])
        self.assertEqual(args.command, "image")
        self.assertEqual(args.image_command, "build")
        self.assertEqual(args.target, "ubuntu:noble")
        self.assertTrue(args.push)

    def test_top_level_image_bake_accepts_component_options(self) -> None:
        parser = build_parser()
        args = parser.parse_args(
            [
                "image",
                "bake",
                "--target",
                "ubuntu:noble",
                "--component",
                "toolboxes",
                "--from-image",
                "ghcr.io/feelpp/feelpp:ubuntu-24.04",
            ]
        )
        self.assertEqual(args.component, "toolboxes")
        self.assertEqual(args.from_image, "ghcr.io/feelpp/feelpp:ubuntu-24.04")

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
            (manifest.parent / "spack.lock").write_text("stale lock\n", encoding="utf-8")
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
            self.assertIn(
                f"spack -e /opt/feelpp/packaging/spack/environments/cpu/openmpi concretize -f",
                dockerfile.read_text(encoding="utf-8"),
            )
            self.assertIn("ARG SPACK_BUILD_JOBS=16", dockerfile.read_text(encoding="utf-8"))
            self.assertIn("ARG SPACK_CONCURRENT_PACKAGES=0", dockerfile.read_text(encoding="utf-8"))
            self.assertIn(
                'spack -e /opt/feelpp/packaging/spack/environments/cpu/openmpi install -j "${SPACK_BUILD_JOBS}";',
                dockerfile.read_text(encoding="utf-8"),
            )
            self.assertIn('"group": {', bake_file.read_text(encoding="utf-8"))
            self.assertIn('"spack-openmpi"', bake_file.read_text(encoding="utf-8"))
            self.assertIn('"ghcr.io/feelpp/feelpp-env:spack-openmpi"', bake_file.read_text(encoding="utf-8"))
            self.assertIn('"packaging_target": "spack:openmpi"', stdout.getvalue())
            self.assertIn('"oci_dist": "spack-openmpi"', stdout.getvalue())
            self.assertIn('"recommended_groups": [', stdout.getvalue())
            self.assertIn('"docker_bake_command": "docker buildx bake -f ', stdout.getvalue())
            self.assertIn(' default"', stdout.getvalue())
            self.assertFalse((context_dir / "packaging" / "spack" / "environments" / "cpu" / "openmpi" / ".spack-env").exists())
            self.assertFalse((context_dir / "packaging" / "spack" / "environments" / "cpu" / "openmpi" / "spack.lock").exists())

    def test_spack_image_bake_uses_manifest_parallelism_and_allows_overrides(self) -> None:
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
            manifest.write_text(
                """
spack:
  config:
    installer: new
    build_jobs: 12
    concurrent_packages: 2
  specs: []
""".strip()
                + "\n",
                encoding="utf-8",
            )

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
                        "--spack-build-jobs",
                        "24",
                        "--spack-concurrent-packages",
                        "3",
                    ]
                )

            context_dir = repo_root / "job" / "images" / "spack-openmpi"
            dockerfile = context_dir / "Dockerfile"
            bake_file = context_dir / "docker-bake.json"
            self.assertEqual(rc, 0)
            self.assertIn("ARG SPACK_BUILD_JOBS=24", dockerfile.read_text(encoding="utf-8"))
            self.assertIn("ARG SPACK_CONCURRENT_PACKAGES=3", dockerfile.read_text(encoding="utf-8"))
            self.assertIn(
                'install -j "${SPACK_BUILD_JOBS}" -p "${SPACK_CONCURRENT_PACKAGES}"',
                dockerfile.read_text(encoding="utf-8"),
            )
            bake_payload = json.loads(bake_file.read_text(encoding="utf-8"))
            self.assertEqual(bake_payload["target"]["spack-openmpi"]["args"]["SPACK_BUILD_JOBS"], "24")
            self.assertEqual(
                bake_payload["target"]["spack-openmpi"]["args"]["SPACK_CONCURRENT_PACKAGES"],
                "3",
            )
            self.assertIn('"spack_build_jobs": 24', stdout.getvalue())
            self.assertIn('"spack_concurrent_packages": 3', stdout.getvalue())

    def test_spack_image_bake_supports_full_build_from_env_image(self) -> None:
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
                        "spack:openmpi",
                        "--component",
                        "full",
                        "--from-image",
                        "ghcr.io/feelpp/feelpp-env:spack-openmpi",
                    ]
                )

            context_dir = Path(tmpdir) / "images" / "spack-openmpi"
            dockerfile = context_dir / "feelpp" / "Dockerfile.full"
            bake_file = context_dir / "docker-bake.json"
            self.assertEqual(rc, 0)
            self.assertTrue(dockerfile.is_file())
            self.assertTrue(bake_file.is_file())
            bake_payload = json.loads(bake_file.read_text(encoding="utf-8"))
            self.assertEqual(bake_payload["group"]["default"]["targets"], ["feelpp-full", "feelpp-full-runtime"])
            self.assertEqual(
                bake_payload["target"]["feelpp-full"]["args"]["FROM_IMAGE"],
                "ghcr.io/feelpp/feelpp-env:spack-openmpi",
            )
            self.assertEqual(
                bake_payload["target"]["feelpp-full"]["args"]["CMAKE_PRESET"],
                "release-clang-spack",
            )
            self.assertEqual(
                bake_payload["target"]["feelpp-full"]["contexts"]["feelpp_source"],
                str(repo_root),
            )
            self.assertIn('"selected_component": "full"', stdout.getvalue())
            self.assertIn('"recommended_groups": [', stdout.getvalue())
            self.assertIn(' full-all', stdout.getvalue())

    def test_top_level_image_targets_lists_repo_owned_images_profile(self) -> None:
        stdout = io.StringIO()
        with mock.patch("sys.stdout", stdout):
            rc = main(["image", "targets", "--repo-root", str(self.repo_root()), "--backend", "apt"])

        self.assertEqual(rc, 0)
        self.assertIn('"target": "ubuntu:noble"', stdout.getvalue())
        self.assertIn('"target": "ubuntu:resolute"', stdout.getvalue())
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
            bake_payload = json.loads(bake_file.read_text(encoding="utf-8"))
            self.assertEqual(bake_payload["group"]["default"]["targets"], ["feelpp-env"])
            self.assertEqual(bake_payload["group"]["env"]["targets"], ["feelpp-env"])
            self.assertIn("feelpp-env", bake_payload["target"])
            self.assertIn("toolboxes-runtime", bake_payload["target"])
            self.assertIn("full-all", bake_payload["group"])
            self.assertEqual(
                bake_payload["target"]["feelpp"]["contexts"]["feelpp_env_image"],
                "target:feelpp-env",
            )
            self.assertEqual(
                bake_payload["target"]["feelpp"]["contexts"]["feelpp_source"],
                str(repo_root),
            )
            self.assertIn("feelpp-all", bake_payload["group"])
            self.assertIn("toolboxes-all", bake_payload["group"])
            self.assertIn("mor-all", bake_payload["group"])
            self.assertIn('"docker_bake_command": "docker buildx bake -f ', stdout.getvalue())
            self.assertIn(' default"', stdout.getvalue())
            self.assertIn('"available_groups": [', stdout.getvalue())
            self.assertIn('"recommended_groups": [', stdout.getvalue())
            self.assertIn('"default_group": "default"', stdout.getvalue())

    def test_top_level_image_bake_supports_external_base_images_for_components(self) -> None:
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
                        "--component",
                        "toolboxes",
                        "--from-image",
                        "ghcr.io/feelpp/feelpp:ubuntu-24.04",
                    ]
                )

            context_dir = Path(tmpdir) / "images" / "ubuntu-noble"
            bake_file = context_dir / "docker-bake.json"
            self.assertEqual(rc, 0)
            bake_payload = json.loads(bake_file.read_text(encoding="utf-8"))
            self.assertEqual(bake_payload["group"]["toolboxes-all"]["targets"], ["toolboxes", "toolboxes-runtime"])
            self.assertEqual(
                bake_payload["target"]["toolboxes"]["args"]["FROM_IMAGE"],
                "ghcr.io/feelpp/feelpp:ubuntu-24.04",
            )
            self.assertEqual(
                bake_payload["target"]["toolboxes"]["contexts"],
                {"feelpp_source": str(repo_root)},
            )
            self.assertEqual(
                bake_payload["target"]["toolboxes"]["args"]["CMAKE_PRESET"],
                "toolboxes",
            )
            self.assertIn('"selected_component": "toolboxes"', stdout.getvalue())

    def test_top_level_image_build_dry_run_uses_load_for_local_builds(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.repo_root()
            stdout = io.StringIO()
            with mock.patch.dict(
                "os.environ",
                {
                    "FEELPP_GITHUB_TOKEN": "",
                    "FEELPP_GIRDER_API_KEY": "",
                    "FEELPP_CKAN_API_KEY": "",
                    "FEELPP_CKAN_URL": "",
                    "FEELPP_CKAN_ORGANIZATION": "",
                },
                clear=False,
            ):
                with mock.patch("sys.stdout", stdout):
                    rc = main(
                        [
                            "image",
                            "build",
                            "--repo-root",
                            str(repo_root),
                            "--job-root",
                            tmpdir,
                            "--target",
                            "ubuntu:noble",
                            "--dry-run",
                        ]
                    )

            bake_file = Path(tmpdir) / "images" / "ubuntu-noble" / "docker-bake.json"
            self.assertEqual(rc, 0)
            self.assertTrue(bake_file.is_file())
            self.assertIn(f"docker buildx bake -f {bake_file} --load default", stdout.getvalue())
            self.assertNotIn("--allow=fs.read=", stdout.getvalue())

    def test_top_level_image_build_dry_run_supports_push_group_override_and_redacts_env(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.repo_root()
            stdout = io.StringIO()
            with mock.patch.dict(
                "os.environ",
                {
                    "FEELPP_GITHUB_TOKEN": "super-secret-token",
                    "FEELPP_CKAN_URL": "https://data.example.invalid",
                },
                clear=False,
            ):
                with mock.patch("sys.stdout", stdout):
                    rc = main(
                        [
                            "image",
                            "build",
                            "--repo-root",
                            str(repo_root),
                            "--job-root",
                            tmpdir,
                            "--target",
                            "ubuntu:noble",
                            "--component",
                            "toolboxes",
                            "--from-image",
                            "ghcr.io/feelpp/feelpp:ubuntu-24.04",
                            "--group",
                            "toolboxes-runtime",
                            "--push",
                            "--dry-run",
                        ]
                    )

            command = stdout.getvalue()
            self.assertEqual(rc, 0)
            self.assertIn(f"--allow=fs.read={repo_root}", command)
            self.assertIn("--push", command)
            self.assertIn("toolboxes-runtime", command)
            self.assertIn("*.args.FEELPP_GITHUB_TOKEN=***", command)
            self.assertIn("*.args.FEELPP_CKAN_URL=***", command)
            self.assertNotIn("super-secret-token", command)
            self.assertNotIn("https://data.example.invalid", command)

    def test_top_level_image_build_passes_fs_read_allow_for_component_builds(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.repo_root()
            with mock.patch("feelpp.pkg.oci.commands.run_checked") as run_checked:
                rc = main(
                    [
                        "image",
                        "build",
                        "--repo-root",
                        str(repo_root),
                        "--job-root",
                        tmpdir,
                        "--target",
                        "ubuntu:noble",
                        "--component",
                        "feelpp",
                        "--from-image",
                        "ghcr.io/feelpp/feelpp-env:ubuntu-24.04",
                    ]
                )

            self.assertEqual(rc, 0)
            run_checked.assert_called_once()
            command = run_checked.call_args.args[0]
            self.assertIn(f"--allow=fs.read={repo_root}", command)

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
