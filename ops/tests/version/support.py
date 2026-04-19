from __future__ import annotations

from pathlib import Path
import json


class VersionRepoMixin:
    def make_repo(self, tmpdir: str) -> Path:
        repo_root = Path(tmpdir) / "repo"
        (repo_root / ".github").mkdir(parents=True)
        (repo_root / "ops").mkdir(parents=True)
        (repo_root / "toolboxes" / "cmake").mkdir(parents=True)
        (repo_root / "mor" / "cmake").mkdir(parents=True)
        (repo_root / "packaging" / "manifest").mkdir(parents=True)
        (repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian").mkdir(parents=True)
        (repo_root / "packaging" / "debian" / "feelpp" / "resolute" / "debian").mkdir(parents=True)
        (repo_root / "packaging" / "debian" / "feelpp-toolboxes" / "noble" / "debian").mkdir(parents=True)
        (repo_root / "packaging" / "debian" / "feelpp-mor" / "noble" / "debian").mkdir(parents=True)
        (repo_root / ".git").mkdir()

        version_template = (
            'set(FEELPP_VERSION_MAJOR "{major}")\n'
            'set(FEELPP_VERSION_MINOR "{minor}")\n'
            'set(FEELPP_VERSION_MICRO "{patch}")\n'
            'set(FEELPP_VERSION_PRERELEASE "{prerelease}")\n'
        )
        for path, values in {
            repo_root / "feelpp.version.cmake": (0, 111, 0, "-preview.13"),
            # Legacy per-component version files may still exist in the tree,
            # but fpp-version should ignore them in favor of the repo root.
            repo_root / "toolboxes" / "cmake" / "feelpp.version.cmake": (0, 108, 0, "-beta.1"),
            repo_root / "mor" / "cmake" / "feelpp.version.cmake": (0, 109, 0, "-beta.1"),
        }.items():
            path.write_text(
                version_template.format(
                    major=values[0],
                    minor=values[1],
                    patch=values[2],
                    prerelease=values[3],
                ),
                encoding="utf-8",
            )

        (repo_root / "codemeta.json").write_text(
            json.dumps(
                {
                    "@context": "https://w3id.org/codemeta/3.0",
                    "@type": "SoftwareSourceCode",
                    "name": "Feel++",
                    "version": "v0.111.0-preview.13",
                    "codeRepository": "https://github.com/feelpp/feelpp",
                    "author": [
                        {"@type": "Person", "givenName": "Christophe", "familyName": "Prud'homme"},
                        {"@type": "Person", "givenName": "Vincent", "familyName": "Chabannes"},
                    ],
                    "contributor": [
                        {"@type": "Person", "givenName": "Thomas", "familyName": "Saigre"},
                    ],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        (repo_root / ".zenodo.json").write_text(
            json.dumps(
                {
                    "title": "Feel++",
                    "version": "v0.111.0-preview.13",
                    "upload_type": "software",
                    "license": "lgpl-3.0",
                    "access_right": "open",
                    "creators": [{"name": "Prud'homme, Christophe"}],
                    "contributors": [{"name": "Saigre, Thomas", "affiliation": "IRMA", "type": "Researcher"}],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        (repo_root / "CITATION.cff").write_text(
            "\n".join(
                [
                    "cff-version: 1.2.0",
                    'message: "If you use this software, please cite it as below."',
                    'title: "Feel++"',
                    "version: v0.111.0-preview.13",
                    "type: software",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

        manifest = "\n".join(
            [
                "version = 1",
                'default_components = ["feelpp", "feelpp-toolboxes", "feelpp-mor"]',
                "",
                "[components.feelpp]",
                'distros = ["noble", "resolute"]',
                "dependencies = []",
                'python_packages = ["python3-feelpp"]',
                "publish = true",
                'package_revision = "2"',
                "",
                '[components."feelpp-toolboxes"]',
                'distros = ["noble"]',
                'dependencies = ["feelpp"]',
                'python_packages = ["python3-feelpp-toolboxes"]',
                "publish = true",
                'package_revision = "4"',
                "",
                '[components."feelpp-mor"]',
                'distros = ["noble"]',
                'dependencies = ["feelpp-toolboxes"]',
                'python_packages = ["python3-feelpp-mor"]',
                "publish = true",
                'package_revision = "5"',
                "",
            ]
        )
        (repo_root / "packaging" / "manifest" / "components.toml").write_text(manifest, encoding="utf-8")
        (repo_root / ".github" / "plan-ci.json").write_text(
            json.dumps(
                {
                    "profiles": {
                        "packaging": {
                            "catalog": {
                                "ubuntu:noble": {"flavor": "ubuntu", "dist": "noble", "version": "24.04"},
                                "ubuntu:resolute": {"flavor": "ubuntu", "dist": "resolute", "version": "26.04"},
                                "debian:trixie": {"flavor": "debian", "dist": "trixie", "version": "13"},
                            }
                        }
                    },
                }
            ),
            encoding="utf-8",
        )
        (repo_root / "ops" / "pyproject.toml").write_text(
            "\n".join(
                [
                    "[project]",
                    'name = "feelpp-ops"',
                    'version = "0.1.0"',
                    'requires-python = ">=3.10"',
                    "",
                    "[tool.feelpp-ops.releaseNotes.publications.hal]",
                    'collections = ["FEEL", "CEMOSIS"]',
                    "rows = 4",
                    'sort = "producedDate_tdate desc"',
                    "",
                ]
            ),
            encoding="utf-8",
        )

        changelog_template = (
            "{source} ({version}) unstable; urgency=medium\n\n"
            "  * Existing entry\n\n"
            " -- Test User <test@example.com>  Mon, 01 Jan 2024 00:00:00 +0000\n"
        )
        changelogs = {
            repo_root / "packaging" / "debian" / "feelpp" / "noble" / "debian" / "changelog": (
                "feelpp",
                "0.111.0~preview.13-1",
            ),
            repo_root / "packaging" / "debian" / "feelpp" / "resolute" / "debian" / "changelog": (
                "feelpp",
                "0.111.0~preview.13-1",
            ),
            repo_root / "packaging" / "debian" / "feelpp-toolboxes" / "noble" / "debian" / "changelog": (
                "feelpp-toolboxes",
                "0.111.0~preview.13-1",
            ),
            repo_root / "packaging" / "debian" / "feelpp-mor" / "noble" / "debian" / "changelog": (
                "feelpp-mor",
                "0.111.0~preview.13-1",
            ),
        }
        for path, (source, version) in changelogs.items():
            path.write_text(changelog_template.format(source=source, version=version), encoding="utf-8")
        return repo_root

