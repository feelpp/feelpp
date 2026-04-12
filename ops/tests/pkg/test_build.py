from __future__ import annotations

from pathlib import Path
import json
import subprocess
import tarfile
import tempfile
import unittest
from unittest import mock
import os

from feelpp.pkg.build import build_component, validate_runtime_linkage
from feelpp.pkg.build.sourcepkg import _prepare_source_tree
from feelpp.pkg.build.outer_prefix import _bootstrap_outer_build_deps, _outer_internal_env
from feelpp.pkg.build.runner import run_pbuilder_build
from feelpp.pkg.config import PackagingContext


class BuildTests(unittest.TestCase):
    def make_context(self, tmpdir: str) -> PackagingContext:
        repo_root = Path(tmpdir) / "repo"
        (repo_root / "packaging" / "pbuilder" / "hooks").mkdir(parents=True)
        (repo_root / "packaging" / "pbuilder" / "pbuilderrc").write_text("", encoding="utf-8")
        auth_script = repo_root / "feelpp" / "tools" / "scripts" / "pkg" / "feelpp_pkg_sudo_auth.sh"
        auth_script.parent.mkdir(parents=True, exist_ok=True)
        auth_script.write_text("#!/bin/sh\n", encoding="utf-8")
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

    def make_packaging_tree(self, context: PackagingContext, component: str = "feelpp") -> Path:
        packaging_dir = context.repo_root / "packaging" / "debian" / component / context.dist / "debian"
        packaging_dir.mkdir(parents=True, exist_ok=True)
        return packaging_dir

    def make_source_archive(self, tmpdir: str, component: str, raw_version: str) -> Path:
        archive_root = Path(tmpdir) / f"{component}-{raw_version}"
        (archive_root / "data").mkdir(parents=True)
        (archive_root / "data" / "payload.txt").write_text("payload\n", encoding="utf-8")
        (archive_root / "data" / "payload-link").symlink_to("payload.txt")

        archive_path = Path(tmpdir) / f"{component}-{raw_version}.tar.gz"
        with tarfile.open(archive_path, "w:gz") as handle:
            handle.add(archive_root, arcname=archive_root.name, recursive=True)
        return archive_path

    def test_prepare_source_tree_flattens_stable_archive_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            self.make_packaging_tree(context)
            archive_path = self.make_source_archive(tmpdir, "feelpp", "0.111.0")
            source_root = context.job_root / "source-packages" / "feelpp"

            with mock.patch("feelpp.pkg.build.sourcepkg.run") as run:
                tree_root, dsc_path, version = _prepare_source_tree(context, "feelpp", archive_path)

            self.assertEqual(version, "0.111.0")
            self.assertEqual(tree_root, source_root / "feelpp-0.111.0")
            self.assertEqual(dsc_path, source_root / "feelpp_0.111.0-1.dsc")
            self.assertTrue((tree_root / "data" / "payload.txt").is_file())
            self.assertTrue((tree_root / "data" / "payload-link").is_symlink())
            self.assertFalse((tree_root / "feelpp-0.111.0").exists())
            run.assert_has_calls(
                [
                    mock.call(
                        [
                            "dch",
                            "-v",
                            "0.111.0-1",
                            "--distribution",
                            "unstable",
                            "-b",
                            "New upstream commits",
                        ],
                        cwd=tree_root,
                    ),
                    mock.call(["dpkg-source", "-b", str(tree_root)], cwd=source_root),
                ]
            )

    def test_prepare_source_tree_flattens_prerelease_archive_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            self.make_packaging_tree(context)
            archive_path = self.make_source_archive(tmpdir, "feelpp", "0.111.0-preview.13")
            source_root = context.job_root / "source-packages" / "feelpp"

            with mock.patch("feelpp.pkg.build.sourcepkg.run") as run:
                tree_root, dsc_path, version = _prepare_source_tree(context, "feelpp", archive_path)

            self.assertEqual(version, "0.111.0~preview.13")
            self.assertEqual(tree_root, source_root / "feelpp-0.111.0~preview.13")
            self.assertEqual(dsc_path, source_root / "feelpp_0.111.0~preview.13-1.dsc")
            self.assertTrue((tree_root / "data" / "payload.txt").is_file())
            self.assertTrue((tree_root / "data" / "payload-link").is_symlink())
            self.assertFalse((tree_root / "feelpp-0.111.0-preview.13").exists())
            run.assert_has_calls(
                [
                    mock.call(
                        [
                            "dch",
                            "-v",
                            "0.111.0~preview.13-1",
                            "--distribution",
                            "unstable",
                            "-b",
                            "New upstream commits",
                        ],
                        cwd=tree_root,
                    ),
                    mock.call(["dpkg-source", "-b", str(tree_root)], cwd=source_root),
                ]
            )

    def test_run_pbuilder_build_passes_mirror_arguments_explicitly(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            dsc_path = context.job_root / "source-packages" / "feelpp" / "feelpp_0.111.0~preview.13-1.dsc"
            dsc_path.parent.mkdir(parents=True, exist_ok=True)
            dsc_path.write_text("", encoding="utf-8")
            result_dir = context.component_results_dir("feelpp")

            mirrorsite = "http://mirror.example/ubuntu/"
            othermirrors = "deb [trusted=yes] file:///tmp/feelpp-pkg-job/local-repo ./"

            with mock.patch("feelpp.pkg.build.runner.pbuilder_mirrorsite", return_value=mirrorsite):
                with mock.patch("feelpp.pkg.build.runner.pbuilder_othermirrors", return_value=othermirrors):
                    with mock.patch("feelpp.pkg.build.runner.run") as run:
                        run_pbuilder_build(context, dsc_path, result_dir)

            run.assert_called_once()
            command = run.call_args.args[0]
            kwargs = run.call_args.kwargs

            self.assertIn("--mirror", command)
            self.assertIn("--othermirror", command)
            self.assertEqual(command[command.index("--mirror") + 1], mirrorsite)
            self.assertEqual(command[command.index("--othermirror") + 1], othermirrors)
            self.assertEqual(command[-1], str(dsc_path))
            self.assertEqual(kwargs["env_overrides"]["MIRRORSITE"], mirrorsite)
            self.assertEqual(kwargs["env_overrides"]["OTHERMIRROR"], othermirrors)
            self.assertEqual(
                kwargs["env_overrides"]["FEELPP_PBUILDER_ALLOW_PUBLIC_FEELPP_FALLBACK"],
                "false",
            )
            self.assertTrue((context.pbuilder_root / "build").is_dir())
            self.assertTrue(result_dir.is_dir())

    def test_run_pbuilder_build_can_enable_public_feelpp_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            dsc_path = context.job_root / "source-packages" / "feelpp-toolboxes" / "feelpp-toolboxes_0.111.0~preview.13-1.dsc"
            dsc_path.parent.mkdir(parents=True, exist_ok=True)
            dsc_path.write_text("", encoding="utf-8")
            result_dir = context.component_results_dir("feelpp-toolboxes")

            with mock.patch("feelpp.pkg.build.runner.pbuilder_mirrorsite", return_value="http://mirror.example/ubuntu/"):
                with mock.patch("feelpp.pkg.build.runner.pbuilder_othermirrors", return_value="deb http://apt.feelpp.org/ubuntu/noble noble latest"):
                    with mock.patch("feelpp.pkg.build.runner.run") as run:
                        run_pbuilder_build(
                            context,
                            dsc_path,
                            result_dir,
                            allow_public_fallback=True,
                        )

            kwargs = run.call_args.kwargs
            self.assertEqual(
                kwargs["env_overrides"]["FEELPP_PBUILDER_ALLOW_PUBLIC_FEELPP_FALLBACK"],
                "true",
            )

    def test_bootstrap_outer_build_deps_uses_local_repo_for_internal_packages(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            context.local_repo_dir.mkdir(parents=True, exist_ok=True)
            (context.local_repo_dir / "Packages").write_text("", encoding="utf-8")
            (context.local_repo_dir / "libfeelpp-toolboxes1-all-dev_0.111.0~preview.13-1_amd64.deb").write_text(
                "fake-deb",
                encoding="utf-8",
            )
            (context.local_repo_dir / "libfeelpp1_0.111.0~preview.13-1_amd64.deb").write_text(
                "fake-deb",
                encoding="utf-8",
            )

            repo_list = Path(tmpdir) / "feelpp-local-repo.list"
            repo_prefs = Path(tmpdir) / "feelpp-local-repo.pref"
            archives_dir = Path(tmpdir) / "archives"
            archives_dir.mkdir(parents=True, exist_ok=True)
            (archives_dir / "libfeelpp-dev_0.111.0~preview.12-1_amd64.deb").write_text(
                "stale",
                encoding="utf-8",
            )
            outer_prefix = context.job_root / "outer-prefix"

            env = {
                "FEELPP_PKG_IN_CONTAINER": "1",
                "FEELPP_PKG_OUTER_LOCAL_REPO_LIST": str(repo_list),
                "FEELPP_PKG_OUTER_LOCAL_REPO_PREFS": str(repo_prefs),
                "FEELPP_PKG_OUTER_APT_ARCHIVES_DIR": str(archives_dir),
            }

            with mock.patch.dict(os.environ, env, clear=False):
                with mock.patch(
                    "feelpp.pkg.build.outer_prefix.collect_component_build_dependencies",
                    return_value=["cmake"],
                ):
                    with mock.patch(
                        "feelpp.pkg.build.outer_prefix.collect_component_internal_build_dependencies",
                        return_value=["libfeelpp-toolboxes1-all-dev"],
                    ):
                        extracted_archives: list[str] = []

                        def fake_subprocess_run(command: list[str], **_: object) -> subprocess.CompletedProcess[str]:
                            self.assertEqual(command[:2], ["dpkg-deb", "-x"])
                            extracted_archives.append(Path(command[2]).name)
                            target = Path(command[-1]) / "usr" / "share" / "feelpp" / "feel" / "cmake" / "modules"
                            target.mkdir(parents=True, exist_ok=True)
                            (target / "Feel++Config.cmake").write_text("", encoding="utf-8")
                            return subprocess.CompletedProcess(command, 0, "", "")

                        with mock.patch("feelpp.pkg.build.outer_prefix.run") as run:
                            with mock.patch("feelpp.pkg.build.outer_prefix.subprocess.run", side_effect=fake_subprocess_run):
                                _bootstrap_outer_build_deps(context, "feelpp-mor")

            self.assertFalse((archives_dir / "libfeelpp-dev_0.111.0~preview.12-1_amd64.deb").exists())
            self.assertEqual(
                repo_list.read_text(encoding="utf-8"),
                f"deb [trusted=yes] file://{context.local_repo_dir} ./\n",
            )
            prefs_text = repo_prefs.read_text(encoding="utf-8")
            self.assertIn("Package: libfeelpp-toolboxes1-all-dev", prefs_text)
            self.assertIn("Pin: origin apt.feelpp.org", prefs_text)
            self.assertTrue((outer_prefix / "usr" / "share" / "feelpp" / "feel" / "cmake" / "modules" / "Feel++Config.cmake").is_file())
            self.assertCountEqual(
                extracted_archives,
                [
                    "libfeelpp-toolboxes1-all-dev_0.111.0~preview.13-1_amd64.deb",
                    "libfeelpp1_0.111.0~preview.13-1_amd64.deb",
                ],
            )

            self.assertEqual(run.call_count, 2)
            update_cmd = run.call_args_list[0].args[0]
            install_cmd = run.call_args_list[1].args[0]
            self.assertEqual(update_cmd[:2], ["apt-get", "update"])
            self.assertIn("cmake", install_cmd)
            self.assertNotIn("libfeelpp-toolboxes1-all-dev", install_cmd)

    def test_outer_internal_env_exports_feelpp_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            prefix = context.job_root / "outer-prefix" / "usr"
            (prefix / "bin").mkdir(parents=True, exist_ok=True)
            (prefix / "lib" / "x86_64-linux-gnu" / "pkgconfig").mkdir(parents=True, exist_ok=True)
            (prefix / "share" / "pkgconfig").mkdir(parents=True, exist_ok=True)
            (prefix / "lib" / "python3" / "dist-packages").mkdir(parents=True, exist_ok=True)

            env = _outer_internal_env(context)

            self.assertEqual(env["FEELPP_DIR"], str(prefix))
            self.assertTrue(env["PATH"].startswith(str(prefix / "bin")))
            self.assertTrue(env["CMAKE_PREFIX_PATH"].startswith(str(prefix)))
            self.assertIn(str(prefix / "lib" / "x86_64-linux-gnu"), env["LD_LIBRARY_PATH"])
            self.assertIn(str(prefix / "lib" / "x86_64-linux-gnu" / "pkgconfig"), env["PKG_CONFIG_PATH"])
            self.assertIn(str(prefix / "lib" / "python3" / "dist-packages"), env["PYTHONPATH"])

    def test_validate_runtime_linkage_rejects_unversioned_feelpp_needed_entry(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            result_dir = Path(tmpdir)
            runtime_deb = result_dir / "libfeelpp1_0.111.0~preview.13-1_amd64.deb"
            runtime_deb.write_text("", encoding="utf-8")

            def fake_subprocess_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
                if command[0] == "dpkg-deb":
                    extract_root = Path(command[-1])
                    target = extract_root / "usr" / "lib" / "x86_64-linux-gnu" / "libfeelpp.so.1"
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(b"\x7fELF")
                    return subprocess.CompletedProcess(command, 0, "", "")
                if command[0] == "readelf":
                    return subprocess.CompletedProcess(
                        command,
                        0,
                        " 0x0000000000000001 (NEEDED)             Shared library: [libfeelpp_specx.so]\n",
                        "",
                    )
                raise AssertionError(f"Unexpected command: {command}")

            with mock.patch("feelpp.pkg.build.elfcheck.subprocess.run", side_effect=fake_subprocess_run):
                with self.assertRaisesRegex(RuntimeError, r"libfeelpp_specx\.so"):
                    validate_runtime_linkage(result_dir)

    def test_validate_runtime_linkage_ignores_dev_packages(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            result_dir = Path(tmpdir)
            dev_deb = result_dir / "libfeelpp-dev_0.111.0~preview.13-1_amd64.deb"
            dev_deb.write_text("", encoding="utf-8")

            with mock.patch("feelpp.pkg.build.elfcheck.subprocess.run") as run_mock:
                validate_runtime_linkage(result_dir)

            run_mock.assert_not_called()

    def test_validate_runtime_linkage_accepts_versioned_feelpp_needed_entry(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            result_dir = Path(tmpdir)
            runtime_deb = result_dir / "feelpp-tools_0.111.0~preview.13-1_amd64.deb"
            runtime_deb.write_text("", encoding="utf-8")

            def fake_subprocess_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
                if command[0] == "dpkg-deb":
                    extract_root = Path(command[-1])
                    target = extract_root / "usr" / "bin" / "feelpp_toolbox_solid"
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(b"\x7fELF")
                    return subprocess.CompletedProcess(command, 0, "", "")
                if command[0] == "readelf":
                    return subprocess.CompletedProcess(
                        command,
                        0,
                        " 0x0000000000000001 (NEEDED)             Shared library: [libfeelpp_fmi4cpp.so.1]\n",
                        "",
                    )
                raise AssertionError(f"Unexpected command: {command}")

            with mock.patch("feelpp.pkg.build.elfcheck.subprocess.run", side_effect=fake_subprocess_run):
                validate_runtime_linkage(result_dir)

    def test_build_component_records_built_component_in_job_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            with mock.patch("feelpp.pkg.build.prepare_runtime_assets"):
                with mock.patch("feelpp.pkg.build.clear_result_dir"):
                    with mock.patch("feelpp.pkg.build.outer_prefix_module._bootstrap_outer_build_deps"):
                        with mock.patch(
                            "feelpp.pkg.build.archive_module._build_source_archive",
                            return_value=context.job_root / "source.tar.gz",
                        ):
                            with mock.patch(
                                "feelpp.pkg.build.sourcepkg_module._prepare_source_tree",
                                return_value=(
                                    context.job_root / "source-packages" / "feelpp",
                                    context.job_root / "source-packages" / "feelpp" / "feelpp_0.111.0~preview.13-1.dsc",
                                    "0.111.0~preview.13",
                                ),
                            ):
                                with mock.patch("feelpp.pkg.build.run_pbuilder_build"):
                                    with mock.patch(
                                        "feelpp.pkg.build.collect_component_internal_build_dependencies",
                                        return_value=[],
                                    ):
                                        with mock.patch("feelpp.pkg.build.validate_runtime_linkage"):
                                            with mock.patch("feelpp.pkg.build.stage_outputs"):
                                                build_component(
                                                    context,
                                                    "feelpp",
                                                skip_pbuilder_prepare=True,
                                            )

            payload = json.loads(context.job_manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["state"], "built")
            self.assertEqual(payload["built_components"], ["feelpp"])


if __name__ == "__main__":
    unittest.main()
