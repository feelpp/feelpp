from __future__ import annotations

import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from feelpp.ops.ci.planner_directives import PlannerDirectiveError, build_planner_message, main


class PlannerDirectiveTests(unittest.TestCase):
    def test_build_message_promotes_spack_only_targets_to_full_mode(self) -> None:
        self.assertEqual(
            build_planner_message(targets="spack:openmpi"),
            "targets=spack:openmpi\nmode=full",
        )

    def test_build_message_normalizes_and_deduplicates_targets(self) -> None:
        self.assertEqual(
            build_planner_message(targets="Spack, openmpi spack-openmpi ubuntu:noble"),
            "targets=spack:openmpi,ubuntu:noble",
        )

    def test_build_message_normalizes_only_skip_and_mode_values(self) -> None:
        self.assertEqual(
            build_planner_message(
                targets="ubuntu:noble",
                only="FeelPP testsuite",
                skip=" MOR ",
                mode="COMPONENTS",
            ),
            "targets=ubuntu:noble\nonly=feelpp,testsuite\nskip=mor\nmode=components",
        )

    def test_build_message_extracts_inline_dispatch_directives_from_targets(self) -> None:
        self.assertEqual(
            build_planner_message(targets="ubuntu:noble,only=testsuite"),
            "targets=ubuntu:noble\nonly=testsuite",
        )

    def test_build_message_extracts_prefixed_inline_dispatch_directives_from_targets(self) -> None:
        self.assertEqual(
            build_planner_message(targets="targets=ubuntu:noble,only=feelpp,testsuite,skip=mor,mode=components"),
            "targets=ubuntu:noble\nonly=feelpp,testsuite\nskip=mor\nmode=components",
        )

    def test_build_message_rejects_component_jobs_for_spack_target(self) -> None:
        with self.assertRaisesRegex(PlannerDirectiveError, "feelpp-full"):
            build_planner_message(targets="spack:openmpi", only="feelpp")

    def test_build_message_rejects_skip_filters_for_spack_target(self) -> None:
        with self.assertRaisesRegex(PlannerDirectiveError, "does not support skip="):
            build_planner_message(targets="spack:openmpi", skip="mor")

    def test_build_message_rejects_non_full_mode_for_spack_only_target(self) -> None:
        with self.assertRaisesRegex(PlannerDirectiveError, "mode=full"):
            build_planner_message(targets="spack:openmpi", mode="components")

    def test_build_message_allows_components_mode_for_mixed_targets(self) -> None:
        self.assertEqual(
            build_planner_message(targets="ubuntu:noble spack:openmpi", mode="components"),
            "targets=ubuntu:noble,spack:openmpi\nmode=components",
        )

    def test_main_writes_message_to_github_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = f"{tmpdir}/github_output.txt"
            rc = main(
                [
                    "--targets",
                    "spack-openmpi",
                    "--github-output",
                    output_path,
                ]
            )

            self.assertEqual(rc, 0)
            contents = Path(output_path).read_text(encoding="utf-8")
            self.assertIn("message<<FEELPP_OUTPUT_", contents)
            self.assertIn("targets=spack:openmpi\nmode=full", contents)

    def test_main_reads_inputs_from_environment(self) -> None:
        stdout = io.StringIO()
        with mock.patch.dict(
            "os.environ",
            {
                "RAW_TARGETS": "ubuntu:noble",
                "RAW_ONLY": "feelpp testsuite",
                "RAW_MODE": "components",
                "GITHUB_OUTPUT": "",
            },
            clear=True,
        ):
            with mock.patch("sys.stdout", stdout):
                rc = main([])

        self.assertEqual(rc, 0)
        self.assertEqual(
            stdout.getvalue(),
            "targets=ubuntu:noble\nonly=feelpp,testsuite\nmode=components\n",
        )

    def test_main_reports_validation_errors_to_stderr(self) -> None:
        stderr = io.StringIO()
        with mock.patch("sys.stderr", stderr):
            rc = main(["--targets", "spack:openmpi", "--only", "feelpp"])

        self.assertEqual(rc, 1)
        self.assertIn("spack:openmpi only supports the feelpp-full job", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
