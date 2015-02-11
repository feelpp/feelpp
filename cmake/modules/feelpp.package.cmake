###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2013-02-08
#
#  Copyright (C) 2013 Université de Strasbourg
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
##
## Archive generation using cpack
##
if (UNIX)
  execute_process(
    COMMAND uname -m
    OUTPUT_VARIABLE FEELPP_SYSTEM_MACHINE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE FEELPP_SYSTEM_MACHINE_error
    RESULT_VARIABLE FEELPP_SYSTEM_MACHINE_result)
endif()
SET(CPACK_PACKAGE_NAME "feelpp")
SET(CPACK_GENERATOR "TGZ")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Feel++")
SET(CPACK_PACKAGE_VENDOR "Christophe Prud'homme")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.md")
SET(CPACK_PACKAGE_VERSION_MAJOR "${FEELPP_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${FEELPP_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${FEELPP_VERSION_MICRO}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "feel")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "CPackSourceConfig.cmake")
SET(CPACK_SYSTEM_NAME "${FEELPP_OS}-${FEELPP_SYSTEM_MACHINE}")

OPTION(FEELPP_ENABLE_CPACK_OPUS "Enable OPUS packaging (if available) in CPack along with Feel++" ON )
SET(CPACK_PACKAGE_NAME "feelpp")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Feel++ (Finite Element method Embedded Library and language in C++)")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${FEELPP_VERSION_MAJOR}.${FEELPP_VERSION_MINOR}.${FEELPP_VERSION_MICRO}${FEELPP_VERSION_PRERELEASE}${FEELPP_VERSION_METADATA}")


SET(CPACK_SOURCE_STRIP_FILES "")
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
set(CPACK_SOURCE_IGNORE_FILES
  "/\\\\.git/;\\\\.gitignore;/\\\\.svn;"
  "/admin/;/Templates/;"
  "/auto/;/ltxpng/;"
  "/TAGS;/#.*;/.*~$;/*.log$;/.cvsignore;/.bzrignore;/work/;/autom4te.cache/"
#  "${PROJECT_SOURCE_DIR}/contrib/eigen"
  "${PROJECT_SOURCE_DIR}/contrib/cln"
  "${PROJECT_SOURCE_DIR}/contrib/boost"
# "${PROJECT_SOURCE_DIR}/contrib/gmsh"
  "${PROJECT_SOURCE_DIR}/contrib/nt2"
  "${PROJECT_SOURCE_DIR}/contrib/ipopt"
  "${PROJECT_SOURCE_DIR}/contrib/cereal"
  "${PROJECT_SOURCE_DIR}/contrib/exodus-5.24"
  "${PROJECT_SOURCE_DIR}/ports"
  "${PROJECT_SOURCE_DIR}/projects"
  "${PROJECT_SOURCE_DIR}/data/medit"
  "${PROJECT_SOURCE_DIR}/data/gmsh/aneurism"
  "${PROJECT_SOURCE_DIR}/data/gmsh/aorta"
  "${PROJECT_SOURCE_DIR}/data/ensoght"
  "${PROJECT_SOURCE_DIR}/data/gmsh/pelvis"
  "${PROJECT_SOURCE_DIR}/applications/"
  "${PROJECT_SOURCE_DIR}/testsuite"
  "${PROJECT_SOURCE_DIR}/benchmarks/"
  "${PROJECT_SOURCE_DIR}/benchmarks/navierstokes/"
  "${PROJECT_SOURCE_DIR}/benchmarks/stokes/"
  "${PROJECT_SOURCE_DIR}/benchmarks/ethiersteinman/"
  "${PROJECT_SOURCE_DIR}/benchmarks/kovasznay/"
  "${PROJECT_SOURCE_DIR}/doc/poster/"
  "${PROJECT_SOURCE_DIR}/doc/api/"
  "${PROJECT_SOURCE_DIR}/doc/manual/pdfs/"
  "${PROJECT_SOURCE_DIR}/doc/manual/pngs/"
  "${PROJECT_SOURCE_DIR}/doc/manual/feel-manual.pdf"
  "${PROJECT_SOURCE_DIR}/doc/manual/manual/feel-manual.pdf"
  "${PROJECT_SOURCE_DIR}/doc/figures/backgrounds/"
  "${PROJECT_SOURCE_DIR}/doc/figures/logos/"
  "${PROJECT_SOURCE_DIR}/doc/manual/feelpp-manual*"
  "${PROJECT_SOURCE_DIR}/doc/manual/*.tex"
  "${PROJECT_SOURCE_DIR}/doc/manual/laplacian/*.msh"
  "${PROJECT_SOURCE_DIR}/doc/manual/heatns/cabin/Mesh/Mesh_NonUniform/*.m*"
  "${PROJECT_SOURCE_DIR}/doc/manual/heatns/cabin/Mesh/Mesh_uniform/*.m*"
  "${PROJECT_SOURCE_DIR}/examples/fluid/"
  "${PROJECT_SOURCE_DIR}/examples/levelset/"
  "${PROJECT_SOURCE_DIR}/examples/pbeq/"
  "${PROJECT_SOURCE_DIR}/research/"
  "/*.tar.gz;*.tar.bz2;*.deb;obj-x86_64-linux-gnu/;"
  "*.eps;*.tag;"
  "*.aux;*.log;*.bbl;*.idx;*.ist;*.out;*.blg;OpusManualBenchmarkEADSUJF.pdf"
  "${PROJECT_SOURCE_DIR}/applications/opus.old"
  "${PROJECT_SOURCE_DIR}/applications/opus/"
   "${PROJECT_SOURCE_DIR}/applications/opus/debian/opus/"
  "${PROJECT_SOURCE_DIR}/applications/opus/debian/source/"
  "${PROJECT_SOURCE_DIR}/applications/opus/doc/"
  "${PROJECT_SOURCE_DIR}/applications/opus/scripts"
  )

#if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/applications/opus/ )
if ( EXISTS ${FEELPP_SOURCE_DIR}/applications/opus/ AND NOT FEELPP_ENABLE_CPACK_OPUS )
#  set(CPACK_SOURCE_IGNORE_FILES "${PROJECT_SOURCE_DIR}/applications/opus/;${CPACK_SOURCE_IGNORE_FILES}" )
endif()

include( CPack )
