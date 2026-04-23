from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.config import (
    DebianPackagingContext,
    PackagingContext,
    WorkspaceContext,
    detect_channel,
    detect_flavor,
)


class ConfigTests(unittest.TestCase):
    def test_detect_channel(self) -> None:
        self.assertEqual(detect_channel("develop"), "latest")
        self.assertEqual(detect_channel("main"), "stable")
        self.assertEqual(detect_channel("feature/pkg"), "latest")

    def test_detect_flavor(self) -> None:
        self.assertEqual(detect_flavor("noble"), "ubuntu")
        self.assertEqual(detect_flavor("resolute"), "ubuntu")
        self.assertEqual(detect_flavor("bookworm"), "debian")

    def test_workspace_context_is_backend_neutral(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir) / "repo"
            context = WorkspaceContext.create(
                repo_root=repo_root,
                branch="develop",
                channel="latest",
                job_id="test-job",
                job_root=repo_root / "job",
            )

            self.assertEqual(context.repo_root, repo_root.resolve())
            self.assertEqual(context.manifest_path, repo_root.resolve() / "packaging" / "manifest" / "components.toml")
            self.assertFalse(hasattr(context, "pbuilder_root"))

    def test_workspace_context_uses_repo_local_default_job_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir) / "repo"
            context = WorkspaceContext.create(
                repo_root=repo_root,
                branch="develop",
                channel="latest",
                job_id="test-job",
            )

            self.assertEqual(
                context.job_root,
                (repo_root / ".cache" / "feelpp-pkg" / "jobs" / "test-job").resolve(),
            )

    def test_packaging_context_alias_points_to_debian_context(self) -> None:
        self.assertIs(PackagingContext, DebianPackagingContext)

    def test_context_uses_runtime_hookdir_under_job_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir) / "repo"
            (repo_root / "packaging" / "pbuilder" / "hooks").mkdir(parents=True)
            (repo_root / "packaging" / "pbuilder" / "pbuilderrc").write_text("", encoding="utf-8")
            job_root = repo_root / "job"
            context = PackagingContext.create(
                repo_root=repo_root,
                dist="noble",
                flavor="ubuntu",
                branch="develop",
                channel="latest",
                job_id="test-job",
                job_root=job_root,
            )

            self.assertEqual(
                context.pbuilder_source_hookdir,
                repo_root / "packaging" / "pbuilder" / "hooks",
            )
            self.assertEqual(
                context.pbuilder_runtime_hookdir,
                job_root / "pbuilder" / "hooks",
            )

    def test_context_uses_distro_scoped_default_pbuilder_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir) / "repo"
            (repo_root / "packaging" / "pbuilder" / "hooks").mkdir(parents=True)
            (repo_root / "packaging" / "pbuilder" / "pbuilderrc").write_text("", encoding="utf-8")
            fake_home = Path(tmpdir) / "home"
            fake_home.mkdir()

            with mock.patch.dict("os.environ", {}, clear=True):
                with mock.patch("pathlib.Path.home", return_value=fake_home):
                    context = PackagingContext.create(
                        repo_root=repo_root,
                        dist="trixie",
                        flavor="debian",
                        branch="develop",
                        channel="latest",
                        job_id="test-job",
                        job_root=repo_root / "job",
                    )

            self.assertEqual(
                context.pbuilder_root,
                fake_home / "pbuilder" / "chroots" / "debian" / "trixie" / "latest",
            )

    def test_context_respects_explicit_pbuilder_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir) / "repo"
            (repo_root / "packaging" / "pbuilder" / "hooks").mkdir(parents=True)
            (repo_root / "packaging" / "pbuilder" / "pbuilderrc").write_text("", encoding="utf-8")
            explicit_root = Path(tmpdir) / "site-pbuilder"

            with mock.patch.dict("os.environ", {"FEELPP_PBUILDER_ROOT": str(explicit_root)}, clear=False):
                context = PackagingContext.create(
                    repo_root=repo_root,
                    dist="noble",
                    flavor="ubuntu",
                    branch="develop",
                    channel="latest",
                    job_id="test-job",
                    job_root=repo_root / "job",
                )

            self.assertEqual(context.pbuilder_root, explicit_root.resolve())


if __name__ == "__main__":
    unittest.main()
