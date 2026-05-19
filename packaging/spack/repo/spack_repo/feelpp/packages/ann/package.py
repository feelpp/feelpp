# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.build_systems.generic import Package


class Ann(Package):
    """Approximate Nearest Neighbor Searching library."""

    homepage = "https://www.cs.umd.edu/~mount/ANN/"
    url = "https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz"

    license("LGPL-2.1-or-later")

    version("1.1.2", sha256="eea03f2e224b66813226d775053316675375dcec45bd263674c052d9324a49a5")

    depends_on("cxx", type="build")

    def install(self, spec, prefix):
        cxx = Executable(spack_cxx)
        mkdirp("spack-build")

        objects = []
        for source in sorted(find("src", "*.cpp")):
            stem = source.replace("/", "_").replace(".cpp", ".o")
            obj = join_path("spack-build", stem)
            cxx("-c", source, "-Iinclude", "-O3", "-fPIC", "-o", obj)
            objects.append(obj)

        cxx("-shared", "-Wl,-soname,libann.so.0", "-o", "libann.so.0.0.0", *objects)
        ar = which("ar")
        ar("crs", "libann.a", *objects)

        mkdirp(prefix.lib)
        install("libann.so.0.0.0", prefix.lib)
        install("libann.a", prefix.lib)
        with working_dir(prefix.lib):
            symlink("libann.so.0.0.0", "libann.so.0")
            symlink("libann.so.0", "libann.so")

        install_tree("include/ANN", prefix.include.ANN)
