from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from feelpp.ops.ci.workflow_plan import compute_workflow_plan, main


def _load_config() -> dict[str, object]:
    return json.loads(
        Path(__file__).resolve().parents[2].joinpath(".github", "plan-ci.json").read_text(encoding="utf-8")
    )


class WorkflowPlanTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = _load_config()

    def test_components_mode_reroutes_spack_to_full(self) -> None:
        outputs = compute_workflow_plan(
            config=self.config,
            mode="components",
            targets=["debian:trixie", "spack:openmpi"],
            enabled_jobs=["feelpp", "testsuite", "toolboxes", "mor"],
        )

        self.assertEqual(json.loads(outputs["component_targets_json"]), ["debian:trixie"])
        self.assertEqual(json.loads(outputs["full_targets_json"]), ["spack:openmpi"])
        self.assertEqual(outputs["run_feelpp"], "true")
        self.assertEqual(outputs["run_testsuite"], "true")
        self.assertEqual(outputs["run_toolboxes"], "true")
        self.assertEqual(outputs["run_mor"], "true")
        self.assertEqual(outputs["run_full"], "true")

    def test_full_mode_uses_full_capable_targets_only(self) -> None:
        outputs = compute_workflow_plan(
            config=self.config,
            mode="full",
            targets=["ubuntu:noble", "spack:openmpi", "debian:trixie"],
            enabled_jobs=["feelpp-full"],
        )

        self.assertEqual(json.loads(outputs["component_targets_json"]), [])
        self.assertEqual(json.loads(outputs["full_targets_json"]), ["ubuntu:noble", "spack:openmpi"])
        self.assertEqual(outputs["run_feelpp"], "false")
        self.assertEqual(outputs["run_toolboxes"], "false")
        self.assertEqual(outputs["run_full"], "true")
        self.assertIn("does not support full builds", json.loads(outputs["warnings_json"])[0])

    def test_component_job_dependencies_are_closed_for_toolboxes(self) -> None:
        outputs = compute_workflow_plan(
            config=self.config,
            mode="components",
            targets=["ubuntu:noble"],
            enabled_jobs=["toolboxes"],
        )

        self.assertEqual(outputs["run_feelpp"], "true")
        self.assertEqual(outputs["run_testsuite"], "false")
        self.assertEqual(outputs["run_toolboxes"], "true")
        self.assertEqual(outputs["run_mor"], "false")
        self.assertEqual(outputs["run_full"], "false")

    def test_component_job_dependencies_are_closed_for_testsuite(self) -> None:
        outputs = compute_workflow_plan(
            config=self.config,
            mode="components",
            targets=["ubuntu:noble"],
            enabled_jobs=["testsuite"],
        )

        self.assertEqual(outputs["run_feelpp"], "true")
        self.assertEqual(outputs["run_testsuite"], "true")
        self.assertEqual(outputs["run_toolboxes"], "false")
        self.assertEqual(outputs["run_mor"], "false")
        self.assertEqual(outputs["run_full"], "false")

    def test_component_job_dependencies_are_closed_for_mor(self) -> None:
        outputs = compute_workflow_plan(
            config=self.config,
            mode="components",
            targets=["ubuntu:noble"],
            enabled_jobs=["mor"],
        )

        self.assertEqual(outputs["run_feelpp"], "true")
        self.assertEqual(outputs["run_toolboxes"], "true")
        self.assertEqual(outputs["run_mor"], "true")

    def test_main_writes_outputs_to_github_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = f"{tmpdir}/github_output.txt"
            rc = main(
                [
                    "--config-path",
                    str(Path(__file__).resolve().parents[2].joinpath(".github", "plan-ci.json")),
                    "--mode",
                    "components",
                    "--targets-json",
                    '["debian:trixie","spack:openmpi"]',
                    "--enabled-jobs-json",
                    '["feelpp","testsuite","toolboxes","mor"]',
                    "--github-output",
                    output_path,
                ]
            )

            self.assertEqual(rc, 0)
            contents = Path(output_path).read_text(encoding="utf-8")
            self.assertIn("component_matrix_json<<FEELPP_OUTPUT_", contents)
            self.assertIn('{"include": [{"target": "debian:trixie"', contents)
            self.assertIn('"target": "spack:openmpi"', contents)


if __name__ == "__main__":
    unittest.main()
