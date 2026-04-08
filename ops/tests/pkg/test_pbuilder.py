from __future__ import annotations

from pathlib import Path
import json
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.config import PackagingContext
from feelpp.pkg.localrepo import stage_outputs
from feelpp.pkg.pbuilder import (
    STAGED_APT_KEYRING_NAME,
    base_tgz_path,
    build_seed_spec,
    pbuilder_othermirrors,
    prepare_base,
    prepare_runtime_assets,
    refresh_feelpp_keyring,
    seed_hash,
    seed_metadata_path,
)


class PbuilderTests(unittest.TestCase):
    def make_context(self, tmpdir: str) -> PackagingContext:
        repo_root = Path(tmpdir) / "repo"
        hooks_dir = repo_root / "packaging" / "pbuilder" / "hooks"
        manifest_dir = repo_root / "packaging" / "manifest"
        control_dir = repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian"
        (hooks_dir / "builddeps").mkdir(parents=True)
        (hooks_dir / "keyrings").mkdir(parents=True)
        manifest_dir.mkdir(parents=True)
        control_dir.mkdir(parents=True)
        (hooks_dir / "G10-feelpp-keyrings").write_text("#!/bin/sh\n", encoding="utf-8")
        (hooks_dir / "H05-feelpp-local-repo-preferences").write_text(
            "#!/bin/sh\n"
            "set -eu\n"
            "cat >/etc/apt/preferences.d/feelpp-local-repo <<'EOF'\n"
            "Package: feelpp-* libfeelpp* python3-feelpp*\n"
            "Pin: origin apt.feelpp.org\n"
            "Pin-Priority: 100\n"
            "EOF\n",
            encoding="utf-8",
        )
        (hooks_dir / "H06-feelpp-local-repo-cache-bust").write_text(
            "#!/bin/sh\n"
            "set -eu\n"
            "find /var/cache/apt/archives -maxdepth 1 -type f \\\n"
            "    \\( -name 'feelpp*.deb' -o -name 'libfeelpp*.deb' -o -name 'python3-feelpp*.deb' \\) \\\n"
            "    -delete\n",
            encoding="utf-8",
        )
        (hooks_dir / "H07-feelpp-local-repo-scrub-stale-files").write_text(
            "#!/bin/sh\n"
            "set -eu\n"
            "rm -rf \\\n"
            "    /usr/include/feelpp \\\n"
            "    /usr/share/feelpp/feel/cmake/modules \\\n"
            "    /usr/include/feelpp/gflags \\\n"
            "    /usr/include/feelpp/glog \\\n"
            "    /usr/include/feelpp/pybind11 \\\n"
            "    /usr/share/feelpp/feelpp_gflags \\\n"
            "    /usr/share/feelpp/gflags \\\n"
            "    /usr/share/feelpp/glog \\\n"
            "    /usr/share/feelpp/pybind11\n"
            "rm -f \\\n"
            "    /usr/share/feelpp/feel/cmake/modules/pybind11*.cmake \\\n"
            "    /usr/share/feelpp/feel/cmake/modules/Findpybind11*.cmake\n",
            encoding="utf-8",
        )
        (hooks_dir / "H08-feelpp-local-repo-apt-update").write_text(
            "#!/bin/sh\n"
            "set -eu\n"
            "local_repo=/tmp/feelpp-pkg-job/local-repo\n"
            "local_repo_list=/etc/apt/sources.list.d/feelpp-local-repo.list\n"
            "public_repo_list=/etc/apt/sources.list.d/feelpp-public.list\n"
            "distribution=${DISTRIBUTION:-noble}\n"
            "flavor=${FLAVOR:-ubuntu}\n"
            "channel=${CHANNEL:-latest}\n"
            "allow_public_fallback=${FEELPP_PBUILDER_ALLOW_PUBLIC_FEELPP_FALLBACK:-false}\n"
            "for source in /etc/apt/sources.list /etc/apt/sources.list.d/*.list /etc/apt/sources.list.d/*.sources; do\n"
            "    [ -f \"${source}\" ] || continue\n"
            "    if ! grep -Fq \"${local_repo}\" \"${source}\"; then\n"
            "        continue\n"
            "    fi\n"
            "    case \"${source}\" in\n"
            "        *.sources)\n"
            "            rm -f \"${source}\"\n"
            "            ;;\n"
            "        *)\n"
            "            tmp=$(mktemp)\n"
            "            grep -Fv \"${local_repo}\" \"${source}\" >\"${tmp}\" || true\n"
            "            if [ -s \"${tmp}\" ]; then\n"
            "                cat \"${tmp}\" >\"${source}\"\n"
            "            else\n"
            "                rm -f \"${source}\"\n"
            "            fi\n"
            "            rm -f \"${tmp}\"\n"
            "            ;;\n"
            "    esac\n"
            "done\n"
            "if [ -f \"${local_repo}/Packages\" ] || [ -f \"${local_repo}/Packages.gz\" ]; then\n"
            "    printf 'deb [trusted=yes] file://%s ./\\n' \"${local_repo}\" >\"${local_repo_list}\"\n"
            "    rm -f \"${public_repo_list}\"\n"
            "elif [ \"${allow_public_fallback}\" = \"1\" ] || [ \"${allow_public_fallback}\" = \"true\" ] || [ \"${allow_public_fallback}\" = \"yes\" ] || [ \"${allow_public_fallback}\" = \"on\" ]; then\n"
            "    printf 'deb http://apt.feelpp.org/%s/%s %s %s\\n' \\\n"
            "        \"${flavor}\" \\\n"
            "        \"${distribution}\" \\\n"
            "        \"${distribution}\" \\\n"
            "        \"${channel}\" \\\n"
            "        >\"${public_repo_list}\"\n"
            "else\n"
            "    rm -f \"${local_repo_list}\"\n"
            "    rm -f \"${public_repo_list}\"\n"
            "fi\n"
            "rm -rf /var/lib/apt/lists/*\n"
            "mkdir -p /var/lib/apt/lists/partial\n"
            "apt-get update\n",
            encoding="utf-8",
        )
        (hooks_dir / "keyrings" / "feelpp.gpg.b64").write_text("ZmFrZQ==\n", encoding="utf-8")
        (repo_root / "packaging" / "pbuilder" / "pbuilderrc").write_text(
            "EXTRAPACKAGES=ca-certificates\n",
            encoding="utf-8",
        )
        (manifest_dir / "components.toml").write_text(
            "\n".join(
                [
                    "version = 1",
                    'default_components = ["feelpp"]',
                    "",
                    "[components.feelpp]",
                    'distros = ["noble"]',
                    "dependencies = []",
                    "python_packages = []",
                    "publish = true",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        (control_dir / "control").write_text(
            "\n".join(
                [
                    "Source: feelpp",
                    "Build-Depends: cmake, debhelper (>= 10), python3:any",
                    "Standards-Version: 3.9.4",
                    "",
                    "Package: libfeelpp-dev",
                    "Architecture: amd64",
                    "Depends: ${misc:Depends}",
                    "Description: dev package",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        pbuilder_root = Path(tmpdir) / "state" / "chroots" / "latest"
        with mock.patch.dict(
            "os.environ",
            {"FEELPP_PBUILDER_ROOT": str(pbuilder_root)},
            clear=False,
        ):
            return PackagingContext.create(
                repo_root=repo_root,
                dist="noble",
                flavor="ubuntu",
                branch="develop",
                channel="latest",
                job_id="test-job",
                job_root=repo_root / "job",
            )

    def test_seed_hash_changes_when_builddeps_change(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            initial_hash = seed_hash(build_seed_spec(context))
            control_file = (
                context.repo_root
                / "packaging"
                / "debian"
                / "feelpp"
                / "noble"
                / "debian"
                / "control"
            )
            control_file.write_text(
                "\n".join(
                    [
                        "Source: feelpp",
                        "Build-Depends: cmake, ninja-build, debhelper (>= 10), python3:any",
                        "Standards-Version: 3.9.4",
                        "",
                        "Package: libfeelpp-dev",
                        "Architecture: amd64",
                        "Depends: ${misc:Depends}",
                        "Description: dev package",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            updated_hash = seed_hash(build_seed_spec(context))

            self.assertNotEqual(initial_hash, updated_hash)

    def test_prepare_runtime_assets_copies_local_repo_preference_hook(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            prepare_runtime_assets(context)

            keyring_hook = context.pbuilder_runtime_hookdir / "G10-feelpp-keyrings"
            self.assertTrue(keyring_hook.is_file())
            self.assertIn("/etc/apt/trusted.gpg.d/feelpp.gpg", keyring_hook.read_text(encoding="utf-8"))
            runtime_hook = context.pbuilder_runtime_hookdir / "H05-feelpp-local-repo-preferences"
            self.assertTrue(runtime_hook.is_file())
            self.assertIn("Pin: origin apt.feelpp.org", runtime_hook.read_text(encoding="utf-8"))
            cache_bust_hook = context.pbuilder_runtime_hookdir / "H06-feelpp-local-repo-cache-bust"
            self.assertTrue(cache_bust_hook.is_file())
            self.assertIn("libfeelpp*.deb", cache_bust_hook.read_text(encoding="utf-8"))
            scrub_hook = context.pbuilder_runtime_hookdir / "H07-feelpp-local-repo-scrub-stale-files"
            self.assertTrue(scrub_hook.is_file())
            self.assertIn("/usr/include/feelpp", scrub_hook.read_text(encoding="utf-8"))
            self.assertIn("/usr/share/feelpp/feel/cmake/modules/pybind11*.cmake", scrub_hook.read_text(encoding="utf-8"))
            apt_refresh_hook = context.pbuilder_runtime_hookdir / "H08-feelpp-local-repo-apt-update"
            self.assertTrue(apt_refresh_hook.is_file())
            self.assertIn("feelpp-local-repo.list", apt_refresh_hook.read_text(encoding="utf-8"))
            self.assertIn("feelpp-public.list", apt_refresh_hook.read_text(encoding="utf-8"))
            self.assertIn("apt.feelpp.org", apt_refresh_hook.read_text(encoding="utf-8"))
            self.assertIn("grep -Fq \"${local_repo}\"", apt_refresh_hook.read_text(encoding="utf-8"))
            self.assertIn("file://%s", apt_refresh_hook.read_text(encoding="utf-8"))
            self.assertIn("apt-get update", apt_refresh_hook.read_text(encoding="utf-8"))

    def test_prepare_runtime_assets_refreshes_feelpp_keyring_from_local_gpg(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            def fake_gpg_export(command: list[str], **_: object) -> subprocess.CompletedProcess[str]:
                output_path = Path(command[command.index("--output") + 1])
                output_path.write_bytes(b"fresh-keyring")
                return subprocess.CompletedProcess(command, 0)

            with mock.patch("feelpp.pkg.pbuilder.assets.subprocess.run", side_effect=fake_gpg_export) as run_mock:
                with mock.patch.dict("os.environ", {"GPG_KEY": "NEWKEY"}, clear=False):
                    prepare_runtime_assets(context)

            keyring_path = context.pbuilder_keyrings_dir / "feelpp.gpg"
            self.assertEqual(b"fresh-keyring", keyring_path.read_bytes())
            run_mock.assert_called_once()
            self.assertIn("NEWKEY", run_mock.call_args.args[0])

    def test_prepare_runtime_assets_keeps_bundled_feelpp_keyring_when_local_export_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            def failing_export(_: list[str], **__: object) -> None:
                raise subprocess.CalledProcessError(2, ["gpg"])

            with mock.patch("feelpp.pkg.pbuilder.assets.run_checked", side_effect=failing_export):
                prepare_runtime_assets(context)

            keyring_path = context.pbuilder_keyrings_dir / "feelpp.gpg"
            self.assertEqual(b"fake", keyring_path.read_bytes())
            keyring_hook = context.pbuilder_runtime_hookdir / "G10-feelpp-keyrings"
            self.assertIn(
                "/etc/apt/trusted.gpg.d/feelpp.gpg",
                keyring_hook.read_text(encoding="utf-8"),
            )

    def test_prepare_runtime_assets_prefers_staged_feelpp_keyring(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            staged_keyring = context.job_root / STAGED_APT_KEYRING_NAME
            staged_keyring.parent.mkdir(parents=True, exist_ok=True)
            staged_keyring.write_bytes(b"staged-keyring")

            with mock.patch("feelpp.pkg.pbuilder.assets.subprocess.run") as run_mock:
                prepare_runtime_assets(context)

            keyring_path = context.pbuilder_keyrings_dir / "feelpp.gpg"
            self.assertEqual(b"staged-keyring", keyring_path.read_bytes())
            keyring_hook = context.pbuilder_runtime_hookdir / "G10-feelpp-keyrings"
            self.assertIn("c3RhZ2VkLWtleXJpbmc=", keyring_hook.read_text(encoding="utf-8"))
            run_mock.assert_not_called()

    def test_prepare_runtime_assets_recreates_missing_runtime_directories(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            context.pbuilder_keyrings_dir.mkdir(parents=True, exist_ok=True)
            context.pbuilder_runtime_hookdir.mkdir(parents=True, exist_ok=True)
            shutil.rmtree(context.pbuilder_keyrings_dir)
            shutil.rmtree(context.pbuilder_runtime_hookdir)

            def fake_gpg_export(command: list[str], **_: object) -> subprocess.CompletedProcess[str]:
                output_path = Path(command[command.index("--output") + 1])
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_bytes(b"fresh-keyring")
                return subprocess.CompletedProcess(command, 0)

            with mock.patch("feelpp.pkg.pbuilder.assets.subprocess.run", side_effect=fake_gpg_export):
                prepare_runtime_assets(context)

            self.assertTrue((context.pbuilder_keyrings_dir / "feelpp.gpg").is_file())
            self.assertTrue((context.pbuilder_runtime_hookdir / "H05-feelpp-local-repo-preferences").is_file())

    def test_refresh_feelpp_keyring_ignores_missing_output_from_gpg(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            target = Path(tmpdir) / "feelpp.gpg"

            with mock.patch(
                "feelpp.pkg.pbuilder.assets.subprocess.run",
                return_value=subprocess.CompletedProcess(["gpg"], 0),
            ):
                self.assertIsNone(refresh_feelpp_keyring(target))

            self.assertFalse(target.exists())

    def test_pbuilder_othermirrors_uses_local_repo_without_public_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            context.local_repo_dir.mkdir(parents=True, exist_ok=True)
            (context.local_repo_dir / "Packages").write_text("", encoding="utf-8")

            mirrors = pbuilder_othermirrors(context)

            self.assertIn(f"deb [trusted=yes] file://{context.local_repo_dir} ./", mirrors)
            self.assertNotIn("apt.feelpp.org", mirrors)

    def test_pbuilder_othermirrors_allows_public_fallback_when_requested(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            context.local_repo_dir.mkdir(parents=True, exist_ok=True)
            (context.local_repo_dir / "Packages").write_text("", encoding="utf-8")

            with mock.patch.dict(
                "os.environ",
                {"FEELPP_PBUILDER_ALLOW_PUBLIC_FEELPP_FALLBACK": "1"},
                clear=False,
            ):
                mirrors = pbuilder_othermirrors(context)

            self.assertIn(f"deb [trusted=yes] file://{context.local_repo_dir} ./", mirrors)
            self.assertNotIn("apt.feelpp.org", mirrors)

    def test_pbuilder_othermirrors_uses_public_fallback_only_without_local_repo(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            with mock.patch.dict(
                "os.environ",
                {"FEELPP_PBUILDER_ALLOW_PUBLIC_FEELPP_FALLBACK": "1"},
                clear=False,
            ):
                mirrors = pbuilder_othermirrors(context)

            self.assertIn("apt.feelpp.org", mirrors)

    def test_prepare_base_reuses_matching_seed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            spec = build_seed_spec(context)
            base_tgz = base_tgz_path(context)
            metadata_path = seed_metadata_path(context)
            base_tgz.parent.mkdir(parents=True, exist_ok=True)
            base_tgz.write_text("placeholder", encoding="utf-8")
            metadata_path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "seed_hash": seed_hash(spec),
                        "spec": spec,
                    }
                ),
                encoding="utf-8",
            )

            with mock.patch("feelpp.pkg.pbuilder.base_tgz_is_valid", return_value=True):
                with mock.patch(
                    "feelpp.pkg.pbuilder.prepare_module._run_prepare_base_shell_with_mirror"
                ) as run_prepare:
                    prepare_base(context)
                    run_prepare.assert_not_called()

    def test_prepare_base_invalidates_stale_seed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            base_tgz = base_tgz_path(context)
            metadata_path = seed_metadata_path(context)
            base_tgz.parent.mkdir(parents=True, exist_ok=True)
            base_tgz.write_text("placeholder", encoding="utf-8")
            metadata_path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "seed_hash": "stale",
                        "spec": {"stale": True},
                    }
                ),
                encoding="utf-8",
            )

            with mock.patch("feelpp.pkg.pbuilder.base_tgz_is_valid", return_value=True):
                with mock.patch(
                    "feelpp.pkg.pbuilder.prepare_module._run_prepare_base_shell_with_mirror"
                ) as run_prepare:
                    prepare_base(context)
                    run_prepare.assert_called_once()
                    self.assertTrue(metadata_path.exists())

    def test_prepare_base_retries_failed_prepare(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            base_tgz = base_tgz_path(context)

            def run_bash_side_effect(*_args, **_kwargs):
                if not hasattr(run_bash_side_effect, "attempt"):
                    run_bash_side_effect.attempt = 0
                run_bash_side_effect.attempt += 1
                if run_bash_side_effect.attempt == 1:
                    raise subprocess.CalledProcessError(1, ["bash", "-lc", "prepare"])
                base_tgz.parent.mkdir(parents=True, exist_ok=True)
                base_tgz.write_bytes(b"seed")

            with mock.patch("feelpp.pkg.pbuilder.base_tgz_is_valid", side_effect=[False, True]):
                with mock.patch(
                    "feelpp.pkg.pbuilder.prepare_module._run_prepare_base_shell_with_mirror",
                    side_effect=run_bash_side_effect,
                ) as run_prepare:
                    prepare_base(context)
                    self.assertEqual(run_prepare.call_count, 2)
                    self.assertTrue(seed_metadata_path(context).exists())

    def test_stage_outputs_replaces_older_same_name_packages(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            first_result = Path(tmpdir) / "result-first"
            second_result = Path(tmpdir) / "result-second"
            first_result.mkdir(parents=True)
            second_result.mkdir(parents=True)

            (first_result / "libfeelpp-dev_0.111.0~preview.12-1_amd64.deb").write_text(
                "old",
                encoding="utf-8",
            )
            with mock.patch("feelpp.pkg.localrepo.subprocess.check_output", return_value=""):
                stage_outputs(context, first_result)
            self.assertTrue(
                (context.local_repo_dir / "libfeelpp-dev_0.111.0~preview.12-1_amd64.deb").is_file()
            )

            (second_result / "libfeelpp-dev_0.111.0~preview.13-1_amd64.deb").write_text(
                "new",
                encoding="utf-8",
            )
            with mock.patch("feelpp.pkg.localrepo.subprocess.check_output", return_value=""):
                stage_outputs(context, second_result)

            self.assertFalse(
                (context.local_repo_dir / "libfeelpp-dev_0.111.0~preview.12-1_amd64.deb").exists()
            )
            self.assertTrue(
                (context.local_repo_dir / "libfeelpp-dev_0.111.0~preview.13-1_amd64.deb").is_file()
            )


if __name__ == "__main__":
    unittest.main()
