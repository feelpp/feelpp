from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.config import PackagingContext
from feelpp.pkg.docker import default_publish_ref, image_ref_tag, publish_docker_image


class DockerTests(unittest.TestCase):
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

    def test_image_ref_tag_uses_last_tag_component(self) -> None:
        self.assertEqual(image_ref_tag("feelpp:noble-preview.13"), "noble-preview.13")
        self.assertEqual(image_ref_tag("ghcr.io/feelpp/feelpp:noble-preview.13"), "noble-preview.13")
        self.assertEqual(image_ref_tag("ghcr.io/feelpp/feelpp"), "latest")

    def test_default_publish_ref_uses_ghcr_defaults(self) -> None:
        self.assertEqual(
            default_publish_ref("feelpp:noble-preview.13"),
            "ghcr.io/feelpp/feelpp:noble-preview.13",
        )

    def test_publish_docker_image_tags_and_pushes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            with mock.patch("feelpp.pkg.docker.run") as run:
                result = publish_docker_image(
                    context,
                    source_ref="feelpp:noble-preview.13",
                    repository="feelpp/custom",
                    dry_run=True,
                )

            self.assertEqual(result["source_ref"], "feelpp:noble-preview.13")
            self.assertEqual(result["target_ref"], "ghcr.io/feelpp/custom:noble-preview.13")
            self.assertEqual(run.call_count, 2)
            self.assertEqual(run.call_args_list[0].args[0], ["docker", "tag", "feelpp:noble-preview.13", "ghcr.io/feelpp/custom:noble-preview.13"])
            self.assertEqual(run.call_args_list[1].args[0], ["docker", "push", "ghcr.io/feelpp/custom:noble-preview.13"])


if __name__ == "__main__":
    unittest.main()
