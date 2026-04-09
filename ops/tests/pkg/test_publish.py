from __future__ import annotations

from pathlib import Path
import os
import subprocess
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.config import PackagingContext
from feelpp.pkg.graph import BuildPlan, ComponentSpec
from feelpp.pkg.publish import publish_cleanup, publish_snapshot
from feelpp.pkg.workspace import write_job_manifest


def _completed(*, returncode: int = 0) -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(args=[], returncode=returncode)


class PublishTests(unittest.TestCase):
    def make_context(self, tmpdir: str) -> PackagingContext:
        repo_root = Path(tmpdir) / "repo"
        repo_root.mkdir(parents=True, exist_ok=True)
        return PackagingContext.create(
            repo_root=repo_root,
            dist="noble",
            flavor="ubuntu",
            branch="develop",
            channel="latest",
            job_id="test-job",
            job_root=repo_root / "job",
        )

    def make_plan(self) -> BuildPlan:
        return BuildPlan(
            dist="noble",
            components=(
                ComponentSpec(
                    name="feelpp",
                    distros=("noble",),
                    dependencies=(),
                    python_packages=("python3-feelpp",),
                    publish=True,
                ),
                ComponentSpec(
                    name="feelpp-toolboxes",
                    distros=("noble",),
                    dependencies=("feelpp",),
                    python_packages=("python3-feelpp-toolboxes",),
                    publish=True,
                ),
            ),
            skipped_components=(),
            publish_enabled=True,
            skip_text="",
        )

    def test_publish_snapshot_rejects_incomplete_recorded_chain(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            write_job_manifest(
                context,
                state="building",
                plan=self.make_plan(),
                built_components=["feelpp"],
            )

            with mock.patch("feelpp.pkg.shell.subprocess.run") as run_mock:
                with self.assertRaisesRegex(RuntimeError, r"feelpp-toolboxes"):
                    publish_snapshot(context)

            run_mock.assert_not_called()

    def test_publish_snapshot_creates_repo_and_publishes_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            write_job_manifest(
                context,
                state="built",
                plan=self.make_plan(),
                built_components=["feelpp", "feelpp-toolboxes"],
            )
            context.artifacts_dir.mkdir(parents=True, exist_ok=True)
            (context.artifacts_dir / "feelpp-tools_1_amd64.deb").write_text("deb", encoding="utf-8")

            with mock.patch.dict(
                os.environ,
                {
                    "FEELPP_PKG_SNAPSHOT_ID": "test-snapshot",
                    "GPG_KEY": "",
                    "GPG_PASSPHRASE": "",
                },
                clear=False,
            ):
                with mock.patch(
                    "feelpp.pkg.shell.subprocess.run",
                    side_effect=[
                        _completed(returncode=1),
                        _completed(),
                        _completed(),
                        _completed(),
                        _completed(returncode=1),
                        _completed(),
                    ],
                ) as run_mock:
                    publish_snapshot(context)

            commands = [call.args[0] for call in run_mock.call_args_list]
            self.assertEqual(
                commands,
                [
                    ["aptly", "repo", "show", "feelpp-noble-latest"],
                    [
                        "aptly",
                        "repo",
                        "create",
                        "-distribution=noble",
                        "-component=latest",
                        "feelpp-noble-latest",
                    ],
                    [
                        "aptly",
                        "repo",
                        "add",
                        "-force-replace",
                        "feelpp-noble-latest",
                        str(context.artifacts_dir),
                    ],
                    [
                        "aptly",
                        "snapshot",
                        "create",
                        "feelpp-noble-latest-snapshot-test-snapshot",
                        "from",
                        "repo",
                        "feelpp-noble-latest",
                    ],
                    [
                        "aptly",
                        "publish",
                        "show",
                        "noble",
                        "s3:apt.feelpp.org:ubuntu/noble",
                    ],
                    [
                        "aptly",
                        "publish",
                        "snapshot",
                        "-force-overwrite",
                        "-distribution=noble",
                        "-component=latest",
                        "feelpp-noble-latest-snapshot-test-snapshot",
                        "s3:apt.feelpp.org:ubuntu/noble",
                    ],
                ],
            )

    def test_publish_snapshot_switches_existing_publish_with_skip_signing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            write_job_manifest(
                context,
                state="built",
                plan=self.make_plan(),
                built_components=["feelpp", "feelpp-toolboxes"],
            )
            context.artifacts_dir.mkdir(parents=True, exist_ok=True)
            (context.artifacts_dir / "feelpp-tools_1_amd64.deb").write_text("deb", encoding="utf-8")

            with mock.patch.dict(
                os.environ,
                {
                    "FEELPP_PKG_SNAPSHOT_ID": "test-snapshot",
                    "FEELPP_APTLY_SKIP_SIGNING": "true",
                    "GPG_KEY": "",
                    "GPG_PASSPHRASE": "",
                },
                clear=False,
            ):
                with mock.patch(
                    "feelpp.pkg.shell.subprocess.run",
                    side_effect=[
                        _completed(),
                        _completed(),
                        _completed(),
                        _completed(),
                        _completed(),
                    ],
                ) as run_mock:
                    publish_snapshot(context)

            commands = [call.args[0] for call in run_mock.call_args_list]
            self.assertEqual(commands[0], ["aptly", "repo", "show", "feelpp-noble-latest"])
            self.assertEqual(
                commands[-1],
                [
                    "aptly",
                    "publish",
                    "switch",
                    "-force-overwrite",
                    "-skip-signing",
                    "-component=latest",
                    "noble",
                    "s3:apt.feelpp.org:ubuntu/noble",
                    "feelpp-noble-latest-snapshot-test-snapshot",
                ],
            )

    def test_publish_snapshot_uses_gpg_key_and_passphrase_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            write_job_manifest(
                context,
                state="built",
                plan=self.make_plan(),
                built_components=["feelpp", "feelpp-toolboxes"],
            )
            context.artifacts_dir.mkdir(parents=True, exist_ok=True)
            (context.artifacts_dir / "feelpp-tools_1_amd64.deb").write_text("deb", encoding="utf-8")

            with mock.patch.dict(
                os.environ,
                {
                    "FEELPP_PKG_SNAPSHOT_ID": "test-snapshot",
                    "GPG_KEY": "ABCDEF0123456789",
                    "GPG_PASSPHRASE": "secret-passphrase",
                },
                clear=False,
            ):
                with mock.patch(
                    "feelpp.pkg.shell.subprocess.run",
                    side_effect=[
                        _completed(returncode=1),
                        _completed(),
                        _completed(),
                        _completed(),
                        _completed(returncode=1),
                        _completed(),
                    ],
                ) as run_mock:
                    publish_snapshot(context)

            publish_command = run_mock.call_args_list[-1].args[0]
            self.assertIn("-gpg-key=ABCDEF0123456789", publish_command)
            passphrase_arg = next(arg for arg in publish_command if arg.startswith("-passphrase-file="))
            passphrase_file = Path(passphrase_arg.split("=", 1)[1])
            self.assertFalse(passphrase_file.exists())

    def test_publish_snapshot_skips_repo_add_when_no_binaries_exist(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            write_job_manifest(
                context,
                state="built",
                plan=self.make_plan(),
                built_components=["feelpp", "feelpp-toolboxes"],
            )
            context.artifacts_dir.mkdir(parents=True, exist_ok=True)

            with mock.patch.dict(
                os.environ,
                {
                    "FEELPP_PKG_SNAPSHOT_ID": "test-snapshot",
                    "GPG_KEY": "",
                    "GPG_PASSPHRASE": "",
                },
                clear=False,
            ):
                with mock.patch(
                    "feelpp.pkg.shell.subprocess.run",
                    side_effect=[
                        _completed(),
                        _completed(),
                        _completed(returncode=1),
                        _completed(),
                    ],
                ) as run_mock:
                    publish_snapshot(context)

            commands = [call.args[0] for call in run_mock.call_args_list]
            self.assertFalse(any(command[0:3] == ["aptly", "repo", "add"] for command in commands))

    def test_publish_cleanup_removes_only_prereleases_with_matching_final_release(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            repo_search_output = "\n".join(
                [
                    "feelpp_0.111.0~preview.12-1_source",
                    "feelpp-data_0.111.0~preview.12-1_all",
                    "feelpp_0.111.0-1_source",
                    "feelpp_0.111.1~preview.1-1_source",
                    "feelpp_0.111.1~rc.1-1_source",
                    "feelpp-mor_0.110.0~beta.2-1_amd64",
                    "feelpp-mor_0.110.0-1_amd64",
                ]
            )

            with mock.patch.dict(
                os.environ,
                {
                    "FEELPP_PKG_SNAPSHOT_ID": "cleanup-snapshot",
                    "GPG_KEY": "",
                    "GPG_PASSPHRASE": "",
                },
                clear=False,
            ):
                with mock.patch(
                    "feelpp.pkg.shell.subprocess.run",
                    side_effect=[
                        _completed(),
                        subprocess.CompletedProcess(
                            args=[],
                            returncode=0,
                            stdout=repo_search_output,
                            stderr="",
                        ),
                        _completed(),
                        _completed(),
                        _completed(returncode=0),
                        _completed(),
                    ],
                ) as run_mock:
                    publish_cleanup(context)

            commands = [call.args[0] for call in run_mock.call_args_list]
            self.assertEqual(commands[0], ["aptly", "repo", "show", "feelpp-noble-latest"])
            self.assertEqual(commands[1], ["aptly", "repo", "search", "feelpp-noble-latest"])
            self.assertEqual(
                commands[2],
                [
                    "aptly",
                    "repo",
                    "remove",
                    "feelpp-noble-latest",
                    "feelpp_0.111.0~preview.12-1_source",
                    "feelpp-data_0.111.0~preview.12-1_all",
                    "feelpp-mor_0.110.0~beta.2-1_amd64",
                ],
            )
            self.assertNotIn("feelpp_0.111.1~preview.1-1_source", commands[2])
            self.assertNotIn("feelpp_0.111.1~rc.1-1_source", commands[2])
            self.assertEqual(
                commands[-1],
                [
                    "aptly",
                    "publish",
                    "switch",
                    "-force-overwrite",
                    "-component=latest",
                    "noble",
                    "s3:apt.feelpp.org:ubuntu/noble",
                    "feelpp-noble-latest-snapshot-cleanup-cleanup-snapshot",
                ],
            )

    def test_publish_cleanup_is_noop_when_no_completed_release_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            repo_search_output = "\n".join(
                [
                    "feelpp_0.111.1~preview.1-1_source",
                    "feelpp_0.111.1~rc.1-1_source",
                ]
            )

            with mock.patch.dict(
                os.environ,
                {
                    "FEELPP_PKG_SNAPSHOT_ID": "cleanup-snapshot",
                    "GPG_KEY": "",
                    "GPG_PASSPHRASE": "",
                },
                clear=False,
            ):
                with mock.patch(
                    "feelpp.pkg.shell.subprocess.run",
                    side_effect=[
                        _completed(),
                        subprocess.CompletedProcess(
                            args=[],
                            returncode=0,
                            stdout=repo_search_output,
                            stderr="",
                        ),
                    ],
                ) as run_mock:
                    publish_cleanup(context)

            commands = [call.args[0] for call in run_mock.call_args_list]
            self.assertEqual(
                commands,
                [
                    ["aptly", "repo", "show", "feelpp-noble-latest"],
                    ["aptly", "repo", "search", "feelpp-noble-latest"],
                ],
            )


if __name__ == "__main__":
    unittest.main()
