from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from feelpp.pkg.config import PackagingContext
from feelpp.pkg.debian import (
    collect_component_build_dependencies,
    collect_seed_build_dependencies,
)


class DebianControlTests(unittest.TestCase):
    def test_collect_seed_build_dependencies_from_control(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir) / "repo"
            manifest_dir = repo_root / "packaging" / "manifest"
            manifest_dir.mkdir(parents=True)
            (manifest_dir / "components.toml").write_text(
                "\n".join(
                    [
                        "version = 1",
                        'default_components = ["feelpp", "feelpp-toolboxes", "feelpp-mor"]',
                        "",
                        "[components.feelpp]",
                        'distros = ["noble"]',
                        "dependencies = []",
                        "python_packages = []",
                        "publish = true",
                        "",
                        '[components."feelpp-toolboxes"]',
                        'distros = ["noble"]',
                        'dependencies = ["feelpp"]',
                        "python_packages = []",
                        "publish = true",
                        "",
                        '[components."feelpp-mor"]',
                        'distros = ["noble"]',
                        'dependencies = ["feelpp-toolboxes"]',
                        "python_packages = []",
                        "publish = true",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            feelpp_control = repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian"
            toolboxes_control = repo_root / "packaging" / "debian" / "feelpp-toolboxes" / "noble" / "debian"
            mor_control = repo_root / "packaging" / "debian" / "feelpp-mor" / "noble" / "debian"
            feelpp_control.mkdir(parents=True)
            toolboxes_control.mkdir(parents=True)
            mor_control.mkdir(parents=True)

            (feelpp_control / "control").write_text(
                "\n".join(
                    [
                        "Source: feelpp",
                        "Build-Depends: quilt, debhelper (>= 10), cmake (>= 3.28), python3:any, libfoo-dev [amd64], bar | baz,",
                        " libzip-dev,",
                        "# comment line",
                        " libboost-all-dev",
                        "Standards-Version: 3.9.4",
                        "",
                        "Package: libfeelpp-dev",
                        "Architecture: amd64",
                        "Depends: ${misc:Depends}",
                        "Description: dev package",
                        "",
                        "Package: feelpp-tools",
                        "Architecture: amd64",
                        "Depends: ${misc:Depends}",
                        "Description: tools package",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            (toolboxes_control / "control").write_text(
                "\n".join(
                    [
                        "Source: feelpp-toolboxes",
                        "Build-Depends: dh-python, libfeelpp-dev, feelpp-tools, python3-dev",
                        "",
                        "Package: libfeelpp-toolboxes1-all-dev",
                        "Architecture: amd64",
                        "Depends: ${misc:Depends}",
                        "Description: toolboxes dev package",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            (mor_control / "control").write_text(
                "\n".join(
                    [
                        "Source: feelpp-mor",
                        "Build-Depends: libfeelpp-toolboxes1-all-dev, python3, dh-python",
                        "",
                        "Package: python3-feelpp-mor",
                        "Architecture: amd64",
                        "Depends: ${misc:Depends}",
                        "Description: mor python package",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            context = PackagingContext.create(
                repo_root=repo_root,
                dist="noble",
                flavor="ubuntu",
                branch="develop",
                channel="latest",
                job_id="test-job",
                job_root=repo_root / "job",
            )

            deps = collect_seed_build_dependencies(context)

            self.assertEqual(
                deps,
                [
                    "quilt",
                    "debhelper",
                    "cmake",
                    "python3",
                    "libfoo-dev",
                    "bar",
                    "libzip-dev",
                    "libboost-all-dev",
                    "dh-python",
                    "python3-dev",
                ],
            )

    def test_resolute_feelpp_build_deps_include_clang(self) -> None:
        repo_root = Path(__file__).resolve().parents[3]
        context = PackagingContext.create(
            repo_root=repo_root,
            dist="resolute",
            flavor="ubuntu",
            branch="develop",
            channel="latest",
            job_id="test-job",
            job_root=repo_root / "build" / "packaging-test-job",
        )

        deps = collect_component_build_dependencies(
            context,
            "feelpp",
            include_prefix_scope=False,
        )

        self.assertIn("clang", deps)
        self.assertNotIn("gcc", deps)


if __name__ == "__main__":
    unittest.main()
