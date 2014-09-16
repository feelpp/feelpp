###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-08-15
#
#  Copyright (C) 2014 Feel++ Consortium
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
# Find the MKL libraries
#
# Options:
#
#   MKL_STATIC        :   use static linking
#   MKL_MULTI_THREADED:   use multi-threading
#   MKL_SDL           :   Single Dynamic Library interface
#
# This module defines the following variables:
#
#   MKL_FOUND            : True if MKL_INCLUDE_DIR are found
#   MKL_INCLUDE_DIR      : where to find mkl.h, etc.
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.


include(FindPackageHandleStandardArgs)

set(INTEL_ROOT "/opt/intel" CACHE PATH "Folder contains intel libs")
set(MKL_ROOT ${INTEL_ROOT}/mkl CACHE PATH "Folder contains MKL")

# Find include dir
find_path(MKL_INCLUDE_DIR mkl.h
    PATHS ${MKL_ROOT}/include)

# Find include directory
#  There is no include folder under linux
if(WIN32)
    find_path(INTEL_INCLUDE_DIR omp.h
        PATHS ${INTEL_ROOT}/include)
    set(MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR} ${INTEL_INCLUDE_DIR})
endif()

# Find libraries

# Handle suffix
set(_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(WIN32)
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES _dll.lib)
    endif()
else()
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    else()
      if ( APPLE )
        set(CMAKE_FIND_LIBRARY_SUFFIXES .dylib)
      else(APPLE)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
      endif()
    endif()
endif()


# MKL is composed by four layers: Interface, Threading, Computational and RTL

if(MKL_SDL)
    find_library(MKL_LIBRARY mkl_rt
        PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/)

    set(MKL_MINIMAL_LIBRARY ${MKL_LIBRARY})
else()
    ######################### Interface layer #######################
    if(WIN32)
        set(MKL_INTERFACE_LIBNAME mkl_intel_c)
    else()
        set(MKL_INTERFACE_LIBNAME mkl_intel)
    endif()

    find_library(MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBNAME}
      PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/)
    message(STATUS "MKL_INTERFACE_LIBRARY: ${MKL_INTERFACE_LIBRARY}")
    ######################## Threading layer ########################
    if(MKL_MULTI_THREADED)
        set(MKL_THREADING_LIBNAME mkl_intel_thread)
    else()
        set(MKL_THREADING_LIBNAME mkl_sequential)
    endif()

    find_library(MKL_THREADING_LIBRARY ${MKL_THREADING_LIBNAME}
      PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/)
    message(STATUS "MKL_THREADING_LIBRARY ${MKL_THREADING_LIBRARY}")

    ####################### Computational layer #####################
    find_library(MKL_CORE_LIBRARY mkl_core
      PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/)
    message(STATUS "MKL_CORE_LIBRARY ${MKL_CORE_LIBRARY}")

    ############################ RTL layer ##########################
    if(WIN32)
        set(MKL_RTL_LIBNAME libiomp5md)
    else()
        set(MKL_RTL_LIBNAME iomp5)
    endif()
    find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME}
      PATHS ${INTEL_ROOT}/lib ${INTEL_ROOT}/lib/intel64 ${INTEL_ROOT}/compiler/lib/intel64/)
    message(STATUS "MKL_RTL_LIBRARY ${MKL_RTL_LIBRARY}")

    # Added check that we have all the components of the MKL library
    if(NOT MKL_INTERFACE_LIBRARY OR NOT MKL_THREADING_LIBRARY OR NOT MKL_CORE_LIBRARY OR NOT MKL_RTL_LIBRARY)
        set(MKL_LIBRARY "")
        set(MKL_MINIMAL_LIBRARY "")
    else()
        set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})
        set(MKL_MINIMAL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})
    endif()
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})



find_package_handle_standard_args(MKL DEFAULT_MSG
    MKL_INCLUDE_DIR MKL_LIBRARY MKL_MINIMAL_LIBRARY)

if(MKL_FOUND)
  set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
  set(MKL_LIBRARIES ${MKL_LIBRARY})
  set(MKL_MINIMAL_LIBRARIES ${MKL_MINIMAL_LIBRARY})
endif()
