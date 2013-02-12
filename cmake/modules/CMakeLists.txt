###  CMakeLists.txt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-12-31
#
#  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#

FILE(GLOB modules "${CMAKE_CURRENT_SOURCE_DIR}/*.cmake")
# set(modules
#   feelpp.extra.astyle.cmake
#   feelpp.extra.warnings.cmake
#   feelpp.macros.cmake
#   feelpp.directive.cmake

#   FindGooglePerfTools.cmake
#   FindGFLAGS.cmake
#   FindGLOG.cmake
#   FindFeel++.cmake
#   FindEigen2.cmake
#   FindMetis.cmake
#   FindOpenTURNS.cmake
#   FindPETSc.cmake
#   FindSLEPc.cmake
#   FindTrilinos.cmake
#   UseDoxygen.cmake
#   FindMadLib.cmake
#   FindOctave.cmake
#   FindPackageMultipass.cmake
#   FindReadline.cmake
#   FindTBB.cmake
#   ResolveCompilerPaths.cmake
#   UseLATEX.cmake )

INSTALL(FILES ${modules} DESTINATION share/feel/cmake/modules COMPONENT Devel)
