# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.packages.gmsh.package import Gmsh as BuiltinGmsh


class Gmsh(BuiltinGmsh):
    """Feel++ overlay for Gmsh package metadata corrections."""

    # The bundled Spack recipe only adds the GUI OpenGL stack for +fltk, but
    # gmsh@4.13.1~fltk still installs libgmsh.so with direct NEEDED entries for
    # these runtime libraries. Keep this in the overlay until upstream metadata
    # matches the produced shared library.
    depends_on("mesa-glu", when="@4.13.1~fltk", type=("build", "link", "run"))
    depends_on("libxcursor", when="@4.13.1~fltk", type=("build", "link", "run"))
    depends_on("libxinerama", when="@4.13.1~fltk", type=("build", "link", "run"))
    depends_on("libxft", when="@4.13.1~fltk", type=("build", "link", "run"))
