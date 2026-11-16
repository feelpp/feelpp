# Copyright Spack Project Developers. See COPYRIGHT file for details.
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
"""ParMmg package with shared-library selection and relocatable CMake exports."""
from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *


class Parmmg(CMakePackage):
    """Parallel mesh adaptation using provisioned MMG, METIS and MPI."""

    homepage = "https://www.mmgtools.org"
    url = "https://github.com/MmgTools/ParMmg/archive/refs/tags/v1.5.0.tar.gz"
    license("LGPL-3.0-or-later")
    version("1.5.0", sha256="0baec7914e49a26bdbb849ab64dcd92147eff79ac02ef3b2599cb05104901a7a")
    variant("shared", default=True, description="Build shared libraries")
    variant("pic", default=True, description="Build position independent code")
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")
    depends_on("mpi")
    depends_on("mmg@5.8: +private_headers")
    depends_on("metis")
    patch("cmake-config.patch")

    def cmake_args(self):
        """Map package variants to upstream's explicit library switches."""
        return [
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("LIBPARMMG_SHARED", "shared"),
            self.define("LIBPARMMG_STATIC", "~shared" in self.spec),
            self.define("DOWNLOAD_MMG", False),
            self.define("DOWNLOAD_METIS", False),
            self.define("MMG_DIR", self.spec["mmg"].prefix),
            self.define("MMG_INCLUDE_DIRS", self.spec["mmg"].headers.directories),
            self.define("METIS_DIR", self.spec["metis"].prefix),
            self.define("USE_VTK", False),
        ]
