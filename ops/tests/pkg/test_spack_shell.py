# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Current-shell activation, failure propagation and human-readable inspection."""
from __future__ import annotations

import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.cli import main


class SpackShellTests(unittest.TestCase):
    """Check shell behavior without provisioning dependencies or using a real Spack."""

    def test_enable_status_disable_and_failed_switch(self) -> None:
        """Activation changes the caller, keeps its Python, and preserves errors."""
        helper = Path(__file__).resolve().parents[2] / "feelpp-spack.sh"
        with tempfile.TemporaryDirectory(prefix="feelpp shell ") as directory:
            root = Path(directory)
            manifests = root / "packaging/spack/environments"
            for name in ("good env", "fail"):
                manifest = manifests / "cpu" / name / "spack.yaml"
                manifest.parent.mkdir(parents=True)
                manifest.write_text("spack:\n  specs: []\n")
            spack = root / "spack/bin/spack"
            spack.parent.mkdir(parents=True)
            spack.write_text(
                f"#!{sys.executable}\n"
                "import sys, shlex\n"
                "if sys.argv[1:3] == ['env', 'activate']:\n"
                "    if sys.argv[4].endswith('/fail'): sys.exit(3)\n"
                "    print('export SPACK_ENV=' + shlex.quote(sys.argv[4]))\n"
                "elif sys.argv[1:3] == ['env', 'deactivate']:\n"
                "    print('unset SPACK_ENV')\n"
                "else: sys.exit(4)\n"
            )
            spack.chmod(0o755)
            environment = {
                **os.environ,
                "FEELPP_OPS_PYTHON": sys.executable,
                "SPACK_ROOT": str(spack.parent.parent),
                "SPACK_USER_CONFIG_PATH": str(root / "config"),
                "SPACK_USER_CACHE_PATH": str(root / "cache"),
                "PATH": str(spack.parent) + os.pathsep + os.environ["PATH"],
            }
            environment.pop("SPACK_ENV", None)
            script = f'''
set -eu
cd() {{ builtin cd "$@" && builtin pwd; }}
pwd() {{ printf 'custom pwd output\n'; }}
source {shlex.quote(str(helper))}
fpp-spack enable "cpu/good env" --repo-root {shlex.quote(str(root))}
test "$SPACK_ENV" = {shlex.quote(str(manifests / 'cpu/good env'))}
fpp-spack status --json
fpp-spack list
if fpp-spack enable cpu/fail >/dev/null 2>&1; then exit 9; fi
test "$SPACK_ENV" = {shlex.quote(str(manifests / 'cpu/good env'))}
if fpp-spack enable missing >/dev/null 2>&1; then exit 10; fi
test "$SPACK_ENV" = {shlex.quote(str(manifests / 'cpu/good env'))}
fpp-spack disable
test -z "${{SPACK_ENV:-}}"
fpp-spack status --json
'''
            result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script],
                                    env=environment, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('"name": "cpu/good env"', result.stdout)
            self.assertIn("* cpu/good env (view missing)", result.stdout)
            self.assertIn('"active": false', result.stdout)

    def test_status_without_active_environment(self) -> None:
        """An inactive shell remains a successful, explicit status result."""
        output = io.StringIO()
        with mock.patch.dict(os.environ, {}, clear=True), mock.patch("sys.stdout", output):
            self.assertEqual(main(["spack", "status", "--json"]), 0)
        self.assertFalse(json.loads(output.getvalue())["active"])

    def test_enable_without_shell_setup_explains_bootstrap(self) -> None:
        """The console command must not pretend to modify its parent shell."""
        output = io.StringIO()
        with mock.patch("sys.stderr", output):
            self.assertEqual(main(["spack", "enable", "cpu/openmpi5-kokkos"]), 1)
        self.assertIn("source", output.getvalue())
        self.assertIn("ops/feelpp-spack.sh", output.getvalue())
