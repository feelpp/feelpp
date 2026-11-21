# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Completion must inspect partial commands without running Spack or activating environments."""
from __future__ import annotations

import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.backends.spack.completion import complete


class SpackCompletionTests(unittest.TestCase):
    """Check command contexts, dynamic manifests and the real Bash integration."""

    def test_commands_and_options_follow_context(self) -> None:
        """Nested commands and flags are offered only in their parser context."""
        self.assertEqual(complete(["ena"]), ["enable"])
        self.assertIn("status", complete([""]))
        self.assertEqual(complete(["env", "sh"]), ["show"])
        self.assertIn("bake", complete(["image", ""]))
        self.assertEqual(complete(["status", "--j"]), ["--json"])
        self.assertEqual(complete(["enable", "--spack-r"]), ["--spack-root"])
        self.assertEqual(complete(["enable", "--repo-root", ""]), [])

    def test_environment_names_aliases_and_option_values(self) -> None:
        """Complete environments with spaces, explicit repo roots and ambiguous aliases."""
        with tempfile.TemporaryDirectory(prefix="feelpp completion ") as directory:
            root = Path(directory)
            for name in ("cpu/openmpi5-kokkos", "cpu/rocm", "hybrid/rocm", "cpu/good env"):
                manifest = root / "packaging/spack/environments" / name / "spack.yaml"
                manifest.parent.mkdir(parents=True)
                manifest.write_text("spack:\n  specs: []\n")
            with mock.patch.dict(os.environ, {"FEELPP_REPO_ROOT": directory}):
                self.assertEqual(complete(["enable", "cpu/openmpi5-k"]), ["cpu/openmpi5-kokkos"])
                self.assertEqual(complete(["activate", "openmpi5-k"]), ["openmpi5-kokkos"])
                self.assertEqual(complete(["env", "show", "cpu/g"]), ["cpu/good env"])
                self.assertEqual(complete(["doctor", "--env", "rocm"]), [])
                self.assertEqual(complete(["doctor", "--env", "hybrid/r"]), ["hybrid/rocm"])
                self.assertEqual(complete(["enable", "--repo-root", directory, "cpu/g"]), ["cpu/good env"])
                self.assertEqual(complete(["enable", f"--repo-root={directory}", "cpu/g"]), ["cpu/good env"])
                self.assertNotIn("cpu/good env", complete(["enable", "cpu/good env", ""]))

    def test_missing_repository_is_quiet(self) -> None:
        """An incomplete or nonexistent repository path produces no environment suggestions."""
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(complete(["enable", "--repo-root", directory, "cpu/"]), [])

    def test_shell_registration_and_pinned_interpreter(self) -> None:
        """Sourcing enables Tab completion, which survives subsequent PATH changes."""
        helper = Path(__file__).resolve().parents[2] / "feelpp-spack.sh"
        script = f'''
set -eu
source {shlex.quote(str(helper))}
complete -p fpp-spack
export PATH=/unavailable
COMP_WORDS=(fpp-spack enable cpu/openmpi5-k)
COMP_CWORD=2
_fpp_spack_complete
printf '%s\\n' "${{COMPREPLY[@]}}"
'''
        result = subprocess.run(
            ["bash", "--noprofile", "--norc", "-c", script],
            env={**os.environ, "FEELPP_OPS_PYTHON": sys.executable},
            text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("-F _fpp_spack_complete fpp-spack", result.stdout)
        self.assertIn("cpu/openmpi5-kokkos", result.stdout.splitlines())
