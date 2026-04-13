from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.apptainer import (
    as_apptainer_source_ref,
    as_oras_ref,
    default_apptainer_binary,
    default_apptainer_publish_ref,
    default_apptainer_tag,
    publish_apptainer_image,
)
from feelpp.pkg.config import PackagingContext


class ApptainerTests(unittest.TestCase):
    def make_context(self, tmpdir: str) -> PackagingContext:
        repo_root = Path(tmpdir) / "repo"
        (repo_root / "packaging" / "pbuilder" / "hooks").mkdir(parents=True)
        (repo_root / "packaging" / "pbuilder" / "pbuilderrc").write_text("", encoding="utf-8")
        job_root = repo_root / "job"
        return PackagingContext.create(
            repo_root=repo_root,
            dist="noble",
            flavor="ubuntu",
            branch="develop",
            channel="latest",
            job_id="test-job",
            job_root=job_root,
        )

    def test_default_apptainer_tag_appends_sif_suffix(self) -> None:
        self.assertEqual(default_apptainer_tag("feelpp:noble-preview.13"), "noble-preview.13_sif")

    def test_default_apptainer_binary_prefers_known_install_path(self) -> None:
        with mock.patch.dict("os.environ", {}, clear=False):
            with mock.patch("feelpp.pkg.apptainer.shutil.which", return_value=None):
                with mock.patch("feelpp.pkg.apptainer.Path.is_file", return_value=True):
                    self.assertEqual(default_apptainer_binary(), "/opt/apptainer/latest/bin/apptainer")

    def test_default_apptainer_publish_ref_uses_ghcr_defaults(self) -> None:
        self.assertEqual(
            default_apptainer_publish_ref("feelpp:noble-preview.13"),
            "ghcr.io/feelpp/feelpp:noble-preview.13_sif",
        )

    def test_as_oras_ref_wraps_plain_refs(self) -> None:
        self.assertEqual(
            as_oras_ref("ghcr.io/feelpp/feelpp:noble-preview.13_sif"),
            "oras://ghcr.io/feelpp/feelpp:noble-preview.13_sif",
        )
        self.assertEqual(
            as_oras_ref("oras://ghcr.io/feelpp/feelpp:noble-preview.13_sif"),
            "oras://ghcr.io/feelpp/feelpp:noble-preview.13_sif",
        )

    def test_as_apptainer_source_ref_uses_docker_transport_for_remote_refs(self) -> None:
        self.assertEqual(
            as_apptainer_source_ref("ghcr.io/feelpp/feelpp:noble-preview.13"),
            "docker://ghcr.io/feelpp/feelpp:noble-preview.13",
        )
        self.assertEqual(
            as_apptainer_source_ref("feelpp:noble-preview.13"),
            "docker-daemon:feelpp:noble-preview.13",
        )

    def test_publish_apptainer_image_builds_and_pushes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            with mock.patch("feelpp.pkg.apptainer.run") as run:
                with mock.patch("feelpp.pkg.apptainer.default_apptainer_binary", return_value="/opt/apptainer/latest/bin/apptainer"):
                    result = publish_apptainer_image(
                        context,
                        source_ref="feelpp:noble-preview.13",
                        repository="feelpp/custom",
                        dry_run=True,
                    )

            self.assertEqual(result["source_ref"], "feelpp:noble-preview.13")
            self.assertEqual(result["resolved_source_ref"], "docker-daemon:feelpp:noble-preview.13")
            self.assertEqual(result["target_ref"], "ghcr.io/feelpp/custom:noble-preview.13_sif")
            self.assertEqual(result["oras_ref"], "oras://ghcr.io/feelpp/custom:noble-preview.13_sif")
            self.assertEqual(result["apptainer_binary"], "/opt/apptainer/latest/bin/apptainer")
            self.assertTrue(result["output_path"].endswith(".sif"))
            self.assertEqual(run.call_count, 2)
            self.assertEqual(
                run.call_args_list[0].args[0],
                [
                    "/opt/apptainer/latest/bin/apptainer",
                    "build",
                    result["output_path"],
                    "docker-daemon:feelpp:noble-preview.13",
                ],
            )
            self.assertEqual(
                run.call_args_list[1].args[0],
                [
                    "/opt/apptainer/latest/bin/apptainer",
                    "push",
                    result["output_path"],
                    "oras://ghcr.io/feelpp/custom:noble-preview.13_sif",
                ],
            )

    def test_publish_apptainer_image_can_build_from_remote_oci_ref(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            with mock.patch("feelpp.pkg.apptainer.run") as run:
                with mock.patch("feelpp.pkg.apptainer.default_apptainer_binary", return_value="/opt/apptainer/latest/bin/apptainer"):
                    result = publish_apptainer_image(
                        context,
                        source_ref="ghcr.io/feelpp/feelpp:noble-preview.13",
                        repository="feelpp/custom",
                        dry_run=True,
                    )

            self.assertEqual(
                result["resolved_source_ref"],
                "docker://ghcr.io/feelpp/feelpp:noble-preview.13",
            )
            self.assertEqual(
                run.call_args_list[0].args[0],
                [
                    "/opt/apptainer/latest/bin/apptainer",
                    "build",
                    result["output_path"],
                    "docker://ghcr.io/feelpp/feelpp:noble-preview.13",
                ],
            )


if __name__ == "__main__":
    unittest.main()
