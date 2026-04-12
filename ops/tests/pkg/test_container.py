from __future__ import annotations

from pathlib import Path
import json
import os
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.config import PackagingContext
from feelpp.pkg.container import (
    CONTAINER_APTLY_CONFIG,
    CONTAINER_GNUPG_HOME,
    CONTAINER_STAGED_GNUPG_HOME,
    DEFAULT_STAGED_APT_KEYRING,
    default_image_for_context,
    publish_uses_aptly,
    publish_uses_signing,
    run_in_docker,
)
from feelpp.pkg.container.runtime import default_state_root


class ContainerTests(unittest.TestCase):
    def make_context(
        self,
        tmpdir: str,
        *,
        dist: str = "noble",
        flavor: str = "ubuntu",
    ) -> PackagingContext:
        repo_root = Path(tmpdir) / "repo"
        (repo_root / "packaging" / "pbuilder" / "hooks").mkdir(parents=True)
        (repo_root / "packaging" / "pbuilder" / "pbuilderrc").write_text("", encoding="utf-8")
        job_root = repo_root / "job"
        return PackagingContext.create(
            repo_root=repo_root,
            dist=dist,
            flavor=flavor,
            branch="develop",
            channel="latest",
            job_id="test-job",
            job_root=job_root,
        )

    def test_publish_uses_signing_for_publish_snapshot(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=False):
            self.assertTrue(publish_uses_signing(["publish", "snapshot", "--dist", "noble"]))
            self.assertTrue(publish_uses_signing(["publish", "cleanup", "--dist", "noble"]))
            self.assertTrue(publish_uses_signing(["build", "chain", "--publish"]))

    def test_publish_uses_aptly_for_publish_snapshot(self) -> None:
        self.assertTrue(publish_uses_aptly(["publish", "snapshot", "--dist", "noble"]))
        self.assertTrue(publish_uses_aptly(["publish", "cleanup", "--dist", "noble"]))
        self.assertTrue(publish_uses_aptly(["build", "chain", "--publish"]))
        self.assertFalse(publish_uses_aptly(["build", "chain"]))

    def test_publish_uses_signing_respects_skip_signing(self) -> None:
        with mock.patch.dict(os.environ, {"FEELPP_APTLY_SKIP_SIGNING": "true"}, clear=False):
            self.assertFalse(publish_uses_signing(["publish", "snapshot", "--dist", "noble"]))

    def test_default_image_for_context_uses_distro_specific_defaults(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            noble = self.make_context(str(Path(tmpdir) / "noble"), dist="noble", flavor="ubuntu")
            trixie = self.make_context(str(Path(tmpdir) / "trixie"), dist="trixie", flavor="debian")
            resolute = self.make_context(str(Path(tmpdir) / "resolute"), dist="resolute", flavor="ubuntu")

            with mock.patch.dict(os.environ, {}, clear=False):
                self.assertEqual(
                    default_image_for_context(noble),
                    "ghcr.io/feelpp/pkg-env:ubuntu-24.04",
                )
                self.assertEqual(
                    default_image_for_context(trixie),
                    "ghcr.io/feelpp/pkg-env:debian-trixie",
                )
                self.assertEqual(
                    default_image_for_context(resolute),
                    "ghcr.io/feelpp/pkg-env:ubuntu-26.04",
                )

    def test_default_image_for_context_respects_override(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            trixie = self.make_context(tmpdir, dist="trixie", flavor="debian")
            with mock.patch.dict(os.environ, {"FEELPP_PKG_DOCKER_IMAGE": "pkg-env:test"}, clear=False):
                self.assertEqual(default_image_for_context(trixie), "pkg-env:test")

    def test_default_state_root_respects_explicit_override(self) -> None:
        with mock.patch.dict(os.environ, {"FEELPP_PKG_DOCKER_STATE_ROOT": "/tmp/feelpp-state"}, clear=True):
            self.assertEqual(default_state_root(), Path("/tmp/feelpp-state").resolve())

    def test_default_state_root_uses_xdg_cache_home(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_home = Path(tmpdir) / "cache"
            with mock.patch("feelpp.pkg.container.runtime._preferred_site_cache_home", return_value=None):
                with mock.patch.dict(os.environ, {"XDG_CACHE_HOME": str(cache_home)}, clear=True):
                    self.assertEqual(default_state_root(), (cache_home / "feelpp-pkg" / "docker").resolve())

    def test_default_state_root_prefers_cemosis_cache_prefix(self) -> None:
        preferred = Path("/nvme0/cemosis/.cache")
        with mock.patch("feelpp.pkg.container.runtime._preferred_site_cache_home", return_value=preferred):
            with mock.patch.dict(os.environ, {}, clear=True):
                self.assertEqual(default_state_root(), (preferred / "feelpp-pkg" / "docker").resolve())

    def test_default_state_root_uses_home_cache_when_xdg_is_unset(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            home = Path(tmpdir) / "home"
            with mock.patch("feelpp.pkg.container.runtime._preferred_site_cache_home", return_value=None):
                with mock.patch.dict(os.environ, {"HOME": str(home)}, clear=True):
                    self.assertEqual(default_state_root(), (home / ".cache" / "feelpp-pkg" / "docker").resolve())

    def test_run_in_docker_stages_host_aptly_config_for_publish(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            aptly_config = Path(tmpdir) / "aptly.conf"
            aptly_config.write_text(
                json.dumps(
                    {
                        "rootDir": "/data/cemosis/aptly",
                        "S3PublishEndpoints": {
                            "apt.feelpp.org": {
                                "region": "eu-west-3",
                                "bucket": "apt.feelpp.org",
                            }
                        },
                    }
                ),
                encoding="utf-8",
            )

            with mock.patch.dict(os.environ, {"FEELPP_APTLY_CONFIG": str(aptly_config)}, clear=False):
                with mock.patch("feelpp.pkg.container.runner.run") as run:
                    run_in_docker(
                        context,
                        argv=["publish", "snapshot", "--dist", "noble"],
                        image="pkg-env:test",
                    )

            command = run.call_args.args[0]
            self.assertIn(f"FEELPP_APTLY_CONFIG={CONTAINER_APTLY_CONFIG}", command)

            staged = json.loads((context.job_root / "aptly.conf").read_text(encoding="utf-8"))
            self.assertEqual(staged["rootDir"], "/srv/aptly")
            self.assertIn("apt.feelpp.org", staged["S3PublishEndpoints"])

    def test_run_in_docker_stages_host_gnupg_for_signed_publish(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            gnupg_home = Path(tmpdir) / "gnupg"
            gnupg_home.mkdir()
            (gnupg_home / "pubring.kbx").write_text("", encoding="utf-8")
            (gnupg_home / "S.gpg-agent").write_text("", encoding="utf-8")

            with mock.patch.dict(os.environ, {"GNUPGHOME": str(gnupg_home)}, clear=False):
                with mock.patch("feelpp.pkg.container.runner.run") as run:
                    run_in_docker(
                        context,
                        argv=["publish", "snapshot", "--dist", "noble"],
                        image="pkg-env:test",
                    )

            command = run.call_args.args[0]
            bootstrap = command[-1]

            self.assertIn(f"export GNUPGHOME={CONTAINER_GNUPG_HOME}", bootstrap)
            self.assertIn(str(CONTAINER_STAGED_GNUPG_HOME), bootstrap)
            self.assertTrue((context.job_root / "gnupg-home" / "pubring.kbx").is_file())
            self.assertFalse((context.job_root / "gnupg-home" / "S.gpg-agent").exists())

    def test_run_in_docker_stages_host_feelpp_apt_key(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            def fake_gpg_export(command: list[str], **_: object) -> object:
                output_path = Path(command[command.index("--output") + 1])
                output_path.write_bytes(b"fresh-public-key")
                return object()

            with mock.patch("feelpp.pkg.container.staging.subprocess.run", side_effect=fake_gpg_export) as run_mock:
                with mock.patch("feelpp.pkg.container.runner.run") as run:
                    run_in_docker(
                        context,
                        argv=["build", "chain", "--dist", "noble"],
                        image="pkg-env:test",
                    )

            self.assertTrue((context.job_root / DEFAULT_STAGED_APT_KEYRING).is_file())
            self.assertEqual(
                b"fresh-public-key",
                (context.job_root / DEFAULT_STAGED_APT_KEYRING).read_bytes(),
            )
            run_mock.assert_called_once()
            run.assert_called_once()

    def test_run_in_docker_bootstrap_installs_debian_archive_keyring_for_debian_dists(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir, dist="trixie", flavor="debian")

            with mock.patch("feelpp.pkg.container.staging.subprocess.run"):
                with mock.patch("feelpp.pkg.container.runner.run") as run:
                    run_in_docker(
                        context,
                        argv=["build", "component", "feelpp", "--dist", "trixie"],
                        image="pkg-env:test",
                    )

            bootstrap = run.call_args.args[0][-1]
            self.assertIn("required_packages='pkgconf arch-test debian-archive-keyring'", bootstrap)

    def test_run_in_docker_skips_host_gnupg_mount_when_signing_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            gnupg_home = Path(tmpdir) / "gnupg"
            gnupg_home.mkdir()

            env = {
                "GNUPGHOME": str(gnupg_home),
                "FEELPP_APTLY_SKIP_SIGNING": "true",
            }
            with mock.patch.dict(os.environ, env, clear=False):
                with mock.patch("feelpp.pkg.container.runner.run") as run:
                    run_in_docker(
                        context,
                        argv=["publish", "snapshot", "--dist", "noble"],
                        image="pkg-env:test",
                    )

            command = run.call_args.args[0]
            self.assertNotIn(str(CONTAINER_STAGED_GNUPG_HOME), command)
            self.assertFalse((context.job_root / "gnupg-home").exists())

    def test_run_in_docker_uses_ops_src_and_maps_repo_stage_result_dir(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            with mock.patch("feelpp.pkg.container.staging.subprocess.run"):
                with mock.patch("feelpp.pkg.container.runner.run") as run:
                    run_in_docker(
                        context,
                        argv=["repo", "stage", "results", "--dist", "noble"],
                        image="pkg-env:test",
                    )

            command = run.call_args.args[0]
            bootstrap = command[-1]
            self.assertIn("export PYTHONPATH=/work/ops/src", bootstrap)
            self.assertIn(
                "python3 -m feelpp.pkg repo stage /work/results --dist noble --repo-root /work --job-root /tmp/feelpp-pkg-job",
                bootstrap,
            )

    def test_run_in_docker_restores_host_ownership_for_workspace_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            with mock.patch("feelpp.pkg.container.staging.subprocess.run"):
                with mock.patch("feelpp.pkg.container.runner.run") as run:
                    run_in_docker(
                        context,
                        argv=["build", "chain", "--dist", "noble"],
                        image="pkg-env:test",
                    )

            bootstrap = run.call_args.args[0][-1]
            self.assertIn("trap cleanup EXIT", bootstrap)
            self.assertIn("for path in /work/build /tmp/feelpp-pkg-job; do", bootstrap)
            self.assertIn(f'chown -R {os.getuid()}:{os.getgid()} "${{path}}" || true', bootstrap)


if __name__ == "__main__":
    unittest.main()
