# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later

from pathlib import Path
import json
import tempfile
import unittest

from feelpp.ops.version.release import ReleaseService
from feelpp.pkg.core.context import WorkspaceContext
from feelpp.pkg.oci.apt import generate_apt_bake
from feelpp.pkg.oci.catalog import load_image_catalog
from feelpp.pkg.oci.common import branch_image_tag
from feelpp.pkg.oci.spack import generate_spack_bake


REPO_ROOT = Path(__file__).resolve().parents[3]
CI_BRANCH = (
    "epic/2383-dynamicstatic-orders-high-order-geometry-fe-refactor-across-"
    "convexes-c20-concepts-kokkos-cpugpu"
)


class OciTagTests(unittest.TestCase):
    def test_existing_short_tags_are_preserved(self) -> None:
        for branch, middle in (("develop", ""), ("", ""), ("feature/test", "-feature-test")):
            for suffix in ("", "-dev", "-full", "-full-dev"):
                with self.subTest(branch=branch, suffix=suffix):
                    self.assertEqual(
                        branch_image_tag("spack-openmpi5", branch=branch, suffix=suffix),
                        f"spack-openmpi5{middle}{suffix}",
                    )

    def test_tag_length_boundary_and_long_branch_collisions(self) -> None:
        for suffix in ("", "-dev", "-full", "-full-dev"):
            budget = 128 - len("spack-openmpi5-") - len(suffix)
            branch = "a" * budget
            self.assertEqual(
                branch_image_tag("spack-openmpi5", branch=branch, suffix=suffix),
                f"spack-openmpi5-{branch}{suffix}",
            )
            tags = [
                branch_image_tag("spack-openmpi5", branch=branch + tail, suffix=suffix)
                for tail in ("x", "y", "/x", "-x")
            ]
            self.assertEqual(len(set(tags)), len(tags))
            for tag in tags:
                self.assertLessEqual(len(tag), 128)
                self.assertRegex(tag, r"^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$")
                self.assertTrue(tag.startswith("spack-openmpi5-"))
                self.assertTrue(tag.endswith(suffix))
            self.assertEqual(
                tags[0], branch_image_tag("spack-openmpi5", branch=branch + "x", suffix=suffix)
            )

    def test_generated_bake_tags_and_release_lookup_for_long_branches(self) -> None:
        catalog = load_image_catalog(REPO_ROOT)
        for branch in (CI_BRANCH, CI_BRANCH + "-extended" * 20):
            for target_name, generator in (
                ("spack:openmpi5", generate_spack_bake),
                ("ubuntu:noble", generate_apt_bake),
            ):
                with self.subTest(branch=branch, target=target_name), tempfile.TemporaryDirectory() as tmpdir:
                    workspace = WorkspaceContext.create(
                        repo_root=REPO_ROOT, branch=branch, job_root=tmpdir
                    )
                    result = generator(
                        workspace=workspace, target=catalog[target_name], component="full",
                        registry="ghcr.io", namespace="feelpp",
                    )
                    bake = json.loads(Path(result["bake_file"]).read_text())
                    for target in bake["target"].values():
                        for ref in target.get("tags", []):
                            tag = ref.rsplit(":", 1)[1]
                            self.assertRegex(tag, r"^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$")
                    self.assertEqual(bake["target"]["feelpp-full"]["args"]["BRANCH"], branch)
                    self.assertTrue(result["image_refs"]["feelpp-full-dev"].endswith("-full-dev"))
                    if target_name.startswith("spack:"):
                        checks = ReleaseService(REPO_ROOT)._spack_container_checks(
                            (target_name,), branch=branch
                        )
                        self.assertEqual(
                            checks[0].candidate_refs, (result["image_refs"]["feelpp-full"],)
                        )
