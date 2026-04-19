from __future__ import annotations

import io
import unittest
from unittest import mock

from feelpp.pkg.cli import build_parser, main


class CliTests(unittest.TestCase):
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

    def test_main_prints_clean_error_for_expected_packaging_failures(self) -> None:
        stderr = io.StringIO()
        with mock.patch("feelpp.pkg.cli_commands.inspect.context_from_args", return_value=object()):
            with mock.patch("feelpp.pkg.cli_commands.inspect.load_plan", side_effect=ValueError("boom")):
                with mock.patch("sys.stderr", stderr):
                    rc = main(["inspect", "plan", "--dist", "trixie"])

        self.assertEqual(rc, 1)
        self.assertEqual(stderr.getvalue().strip(), "error: boom")


if __name__ == "__main__":
    unittest.main()
