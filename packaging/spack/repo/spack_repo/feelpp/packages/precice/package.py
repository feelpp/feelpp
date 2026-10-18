# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import version
from spack_repo.builtin.packages.precice.package import Precice as BuiltinPrecice


class Precice(BuiltinPrecice):
    """Pin the coupling release used by Feel++ without duplicating the builtin recipe."""

    version("3.4.1", sha256="ef4713c938a1b2000d0b071175e1b45f9ec55c7aec4bbe7b65c3992edcc74ac7")

    def cmake_args(self):
        """Install the version probe without building the upstream test suite."""
        return super().cmake_args() + [self.define("BUILD_TESTING", False),
                                     self.define("PRECICE_BUILD_TOOLS", True)]
