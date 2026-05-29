from __future__ import annotations

import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from feelpp.ops.dev.cli import build_parser, main
from feelpp.ops.dev.devcontainer import (
    default_targets,
    load_profiles,
    render_files,
)


PLAN_CI = """
{
  "profiles": {
    "images": {
      "defaults": {
        "targets": ["ubuntu:noble"]
      },
      "catalog": {
        "ubuntu:noble": {
          "flavor": "ubuntu",
          "dist": "noble",
          "version": "24.04",
          "docker": "true",
          "image_backend": "apt",
          "image_strategy": "environment",
          "base_image": "ubuntu:24.04",
          "oci_dist": "ubuntu-24.04"
        },
        "debian:trixie": {
          "flavor": "debian",
          "dist": "trixie",
          "version": "13",
          "docker": "true",
          "image_backend": "apt",
          "image_strategy": "environment",
          "base_image": "debian:13",
          "oci_dist": "debian-13"
        },
        "spack:openmpi": {
          "flavor": "spack",
          "dist": "openmpi",
          "version": "latest",
          "docker": "true",
          "image_backend": "spack",
          "image_strategy": "environment",
          "base_image": "ubuntu:24.04",
          "oci_dist": "spack-openmpi",
          "cmake_full_preset": "release-clang-spack",
          "spack_environment": "cpu/openmpi"
        }
      },
      "groups": {
        "env-current": ["ubuntu:noble", "debian:trixie", "spack:openmpi"]
      }
    }
  }
}
""".strip()


def make_repo(tmpdir: str) -> Path:
    repo_root = Path(tmpdir)
    plan_path = repo_root / ".github" / "plan-ci.json"
    plan_path.parent.mkdir(parents=True)
    plan_path.write_text(PLAN_CI + "\n", encoding="utf-8")
    return repo_root


class DevcontainerTests(unittest.TestCase):
    def test_parser_exposes_devcontainer_commands(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["devcontainer", "generate", "--all"])
        self.assertEqual(args.command, "devcontainer")
        self.assertEqual(args.devcontainer_command, "generate")
        self.assertTrue(args.all)

    def test_default_targets_come_from_image_group(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = make_repo(tmpdir)
            self.assertEqual(
                default_targets(repo_root),
                ["ubuntu:noble", "debian:trixie", "spack:openmpi"],
            )

    def test_profiles_derive_names_images_and_cmake_presets(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = make_repo(tmpdir)
            profiles = load_profiles(repo_root)
            by_id = {profile.id: profile for profile in profiles}

            self.assertEqual(
                by_id["debian-trixie"].image,
                "ghcr.io/feelpp/feelpp-env:debian-13",
            )
            self.assertEqual(
                by_id["debian-trixie"].cmake_preset,
                "container-debian-trixie",
            )
            self.assertEqual(by_id["debian-trixie"].inherits_preset, "default")
            self.assertEqual(
                by_id["spack-openmpi"].inherits_preset,
                "release-clang-spack",
            )

    def test_render_files_writes_profiles_root_default_and_cmake_presets(self) -> None:
        with tempfile.TemporaryDirectory(prefix="feelpp-") as tmpdir:
            repo_root = make_repo(tmpdir)
            files = render_files(repo_root)
            root_devcontainer = json.loads(files[Path(".devcontainer/devcontainer.json")])
            trixie_devcontainer = json.loads(
                files[Path(".devcontainer/debian-trixie/devcontainer.json")]
            )
            presets = json.loads(files[Path("CMakeUserPresets.json")])

            self.assertEqual(root_devcontainer, trixie_devcontainer)
            self.assertEqual(
                root_devcontainer["image"],
                "ghcr.io/feelpp/feelpp-env:debian-13",
            )
            self.assertEqual(
                root_devcontainer["customizations"]["vscode"]["settings"]["cmake.configurePreset"],
                "container-debian-trixie",
            )
            preset_names = {preset["name"] for preset in presets["configurePresets"]}
            self.assertIn("container-debian-trixie", preset_names)
            self.assertIn("container-ubuntu-noble", preset_names)
            self.assertIn("container-spack-openmpi", preset_names)

    def test_cli_generate_and_validate(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = make_repo(tmpdir)
            stdout = io.StringIO()
            with mock.patch("sys.stdout", stdout):
                rc = main(["--repo-root", str(repo_root), "devcontainer", "generate"])

            self.assertEqual(rc, 0)
            self.assertTrue(
                (repo_root / ".devcontainer" / "debian-trixie" / "devcontainer.json").is_file()
            )
            self.assertTrue(
                (repo_root / ".devcontainer" / "ubuntu-noble" / "devcontainer.json").is_file()
            )
            self.assertTrue(
                (repo_root / ".devcontainer" / "spack-openmpi" / "devcontainer.json").is_file()
            )
            self.assertTrue((repo_root / "CMakeUserPresets.json").is_file())
            self.assertIn(".devcontainer/debian-trixie/devcontainer.json", stdout.getvalue())

            with mock.patch("sys.stdout", io.StringIO()):
                validate_rc = main(["--repo-root", str(repo_root), "devcontainer", "validate"])

            self.assertEqual(validate_rc, 0)

    def test_cli_check_reports_drift_without_writing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = make_repo(tmpdir)
            stdout = io.StringIO()
            with mock.patch("sys.stdout", stdout):
                rc = main(["--repo-root", str(repo_root), "devcontainer", "generate", "--check"])

            self.assertEqual(rc, 1)
            self.assertFalse((repo_root / "CMakeUserPresets.json").exists())
            self.assertIn('"drifted"', stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
