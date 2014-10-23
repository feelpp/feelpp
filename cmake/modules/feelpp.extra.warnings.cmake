# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-03-15
#
#  Copyright (C) 2012 Université Joseph Fourier
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 3.0 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#

################################################################################
# Setup some compliler options
################################################################################
include(CheckCXXCompilerFlag)

################################################################################
# Add Extra warning level
################################################################################
option(FEELPP_EXTRA_WARNINGS "Enable/disable extra warnings" OFF)

if(FEELPP_EXTRA_WARNINGS)
  message(STATUS "[feelpp] extra warnings enabled")
endif()

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(HAS_GCC_WALL 1)
  set(HAS_MSVC_W4 0)
elseif(MSVC)
  set(HAS_GCC_WALL 0)
  set(HAS_GCC_WEXTRA 0)
  set(HAS_MSVC_W4 1)
endif()

################################################################################
# Check for GCC style
################################################################################
if(FEELPP_EXTRA_WARNINGS)
  if(NOT DEFINED HAS_GCC_WEXTRA)
    check_cxx_compiler_flag("-Wextra" HAS_GCC_WEXTRA)
  endif()
  if(HAS_GCC_WEXTRA)
    set(HAS_GCC_WALL 1)
  endif()
endif()

if(NOT DEFINED HAS_GCC_WALL)
  check_cxx_compiler_flag("-Wall" HAS_GCC_WALL)
endif()

if(HAS_GCC_WALL)
  set(FEELPP_FLAGS "${FEELPP_FLAGS} -Wall -Wno-unused -Wno-sign-compare ")
endif()
if(FEELPP_EXTRA_WARNINGS)
  if(HAS_GCC_WEXTRA)
    set(FEELPP_FLAGS "${FEELPP_FLAGS} -Wextra")
  endif()
endif()


if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
      CHECK_CXX_COMPILER_FLAG( "-Wno-deprecated-register" HAS_NO_DEPRECATED_REGISTER )
   endif()
else()
   CHECK_CXX_COMPILER_FLAG( "-Wno-deprecated-register" HAS_NO_DEPRECATED_REGISTER )
endif()

if (  HAS_NO_DEPRECATED_REGISTER )
  set( FEELPP_FLAGS "${FEELPP_FLAGS} -Wno-deprecated-register")
endif()

################################################################################
# Check for MSVC style
################################################################################
if(FEELPP_EXTRA_WARNINGS)
  if(NOT DEFINED HAS_MSVC_W4)
    check_cxx_compiler_flag("/W4" HAS_MSVC_W4)
  endif()

  if(HAS_MSVC_W4)
    set(FEELPP_FLAGS "${FEELPP_FLAGS} /W4")
  endif()
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
   set(FEELPP_FLAGS "${FEELPP_FLAGS} -fmacro-backtrace-limit=0" )
   set(FEELPP_FLAGS "${FEELPP_FLAGS} -ftemplate-backtrace-limit=0" )
endif()
