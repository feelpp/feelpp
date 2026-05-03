from __future__ import annotations

import subprocess
import sys
import unittest

from feelpp.pkg.shell import run_capture, run_probe


class ShellTests(unittest.TestCase):
    def test_run_capture_applies_env_overrides_without_explicit_env(self) -> None:
        output = run_capture(
            [
                sys.executable,
                "-c",
                "import os; print(os.environ['FEELPP_TEST_OVERRIDE'])",
            ],
            env_overrides={"FEELPP_TEST_OVERRIDE": "from-override"},
        )

        self.assertEqual(output.strip(), "from-override")

    def test_run_capture_env_overrides_replace_explicit_env(self) -> None:
        output = run_capture(
            [
                sys.executable,
                "-c",
                "import os; print(os.environ['FEELPP_TEST_OVERRIDE'])",
            ],
            env={"FEELPP_TEST_OVERRIDE": "from-base"},
            env_overrides={"FEELPP_TEST_OVERRIDE": "from-override"},
        )

        self.assertEqual(output.strip(), "from-override")

    def test_run_probe_uses_env_overrides(self) -> None:
        succeeded = run_probe(
            [
                sys.executable,
                "-c",
                "import os, sys; sys.exit(0 if os.environ.get('FEELPP_TEST_FLAG') == '1' else 7)",
            ],
            env_overrides={"FEELPP_TEST_FLAG": "1"},
        )

        self.assertTrue(succeeded)

    def test_run_capture_raises_with_merged_env(self) -> None:
        with self.assertRaises(subprocess.CalledProcessError):
            run_capture(
                [
                    sys.executable,
                    "-c",
                    "import os, sys; sys.exit(0 if os.environ.get('FEELPP_TEST_FLAG') == '1' else 9)",
                ],
                env_overrides={"FEELPP_TEST_FLAG": "0"},
            )


if __name__ == "__main__":
    unittest.main()
