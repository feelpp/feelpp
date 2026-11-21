# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
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
SET( SRC
  exprbase.cpp
  expr.cpp
  ginac.cpp
  ginacbase.cpp
  detail/ginacbuildlibrary.cpp
  detail/ginacjit.cpp
  detail/ginacfunctions.cpp
  bilinearformbase.cpp
  symbolsexpr.cpp
  )

set(FEELPP_JIT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/detail")
set(FEELPP_EXCOMPILER_COMMAND "${CMAKE_CURRENT_BINARY_DIR}/feelpp_excompiler")
set(_USE_STANDARD_C_HEADERS_ "-DUSE_STANDARD_HEADERS_IN_GINAC_EXCOMPILER")
configure_file(detail/feelpp_excompiler.in "${FEELPP_EXCOMPILER_COMMAND}" @ONLY)
file(CHMOD "${FEELPP_EXCOMPILER_COMMAND}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
configure_file(detail/ginacjitconfig.hpp.in "${CMAKE_CURRENT_BINARY_DIR}/ginacjitconfig.hpp" @ONLY)

# Installed compilation uses only installed Feel++ headers.
set(FEELPP_JIT_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include/feelpp/feel/feelvf/detail")
configure_file(detail/feelpp_excompiler.in "${CMAKE_CURRENT_BINARY_DIR}/install/feelpp_excompiler" @ONLY)
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/install/feelpp_excompiler" DESTINATION bin COMPONENT Bin)

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp detail/*.hpp)
feel_append_src(feelvf "${SRC};${HDR}")
