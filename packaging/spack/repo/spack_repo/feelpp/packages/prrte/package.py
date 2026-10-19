# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.prrte.package import Prrte as BuiltinPrrte

from spack.package import depends_on, join_path, patch


class Prrte(BuiltinPrrte):
    """Restore binding warnings missing from the PRRTE 4.1.0 release.

    A denied memory-binding request in an ordinary Docker container otherwise
    makes the launcher's warning formatter dereference a null message.
    """

    patch("restore-binding-warnings.patch", when="@4.1.0")
    depends_on("python@3.7:", type="build", when="@4.1.0")

    def autoreconf(self, spec, prefix):
        """Regenerate embedded messages after patching the release help files.

        The tarball ships a prebuilt C table whose make rule does not depend
        on the help text files. Generate without purging: the 4.1.0 citation
        scanner does not recognize the launcher's warning helper calls.
        """
        super().autoreconf(spec, prefix)
        if spec.satisfies("@4.1.0"):
            source = self.stage.source_path
            spec["python"].command(
                join_path(source, "src", "util", "prte-convert-help.py"),
                "--root", source,
                "--out", join_path(source, "src", "util", "prte_show_help_content.c"),
            )
