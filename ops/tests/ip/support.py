from __future__ import annotations

from pathlib import Path
from typing import Any

from feelpp.ops.ip.yamlio import write_yaml


def minimal_public_metadata() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "software": {
            "name": "Feel++",
            "short_name": "feelpp",
            "description": "Public test metadata.",
            "homepage": "https://docs.feelpp.org",
        },
        "repository": {
            "hosting": "GitHub",
            "owner": "feelpp",
            "name": "feelpp",
            "url": "https://github.com/feelpp/feelpp",
            "default_branch": "develop",
        },
        "languages": [{"name": "C++"}, {"name": "Python"}],
        "build_tools": [{"name": "CMake"}, {"name": "Setuptools"}],
        "public_license": {
            "expressions": ["LGPL-3.0-or-later", "GPL-3.0-or-later"],
            "files": [{"path": "LICENSE"}],
        },
        "architecture": {
            "overview": "Public architecture overview.",
            "components": [
                {
                    "name": "C++ core",
                    "paths": ["feelpp/feel"],
                    "description": "Shared C++ library.",
                }
            ],
        },
        "public_funding": {"grants": []},
        "public_outputs": {
            "citation": {"path": "CITATION.cff"},
            "codemeta": {"path": "codemeta.json"},
        },
        "metrics": {
            "approximate_lines_of_code": None,
            "approximate_bytes": None,
            "counted_files": None,
            "method": "tracked text source files",
        },
        "private_boundary": {
            "private_dossier_repository": "cemosis/software-ip",
            "excluded_information": [
                "private inventor records",
                "HR data",
                "signature material",
                "private ownership allocations",
            ],
        },
    }


def write_public_metadata(repo_root: Path, payload: dict[str, Any] | None = None) -> Path:
    metadata_path = repo_root / "metadata" / "software.public.yml"
    write_yaml(metadata_path, payload or minimal_public_metadata())
    return metadata_path
