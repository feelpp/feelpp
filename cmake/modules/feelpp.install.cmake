###  coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2013-02-08
#
#  Copyright (C) 2013-2015 Feel++ Consortium
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
#-----------------------------------------------------------------------
# Values whose defaults are relative to DATAROOTDIR.
# Store empty values in the cache and store the defaults in local
# variables if the cache values are not set explicitly. This auto-updates
# the defaults as DATAROOTDIR changes.
#
if(NOT CMAKE_INSTALL_DATADIR)
  set(CMAKE_INSTALL_DATADIR "" CACHE PATH "read-only architecture-independent data (DATAROOTDIR)")
  set(CMAKE_INSTALL_DATADIR "${CMAKE_INSTALL_DATAROOTDIR}")
endif()

if(NOT CMAKE_INSTALL_MANDIR)
  set(CMAKE_INSTALL_MANDIR "" CACHE PATH "man documentation (DATAROOTDIR/man)")
  set(CMAKE_INSTALL_MANDIR "${CMAKE_INSTALL_DATAROOTDIR}/man")
endif()

if(NOT CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR "" CACHE PATH "documentation root (DATAROOTDIR/doc/PROJECT_NAME)")
  set(CMAKE_INSTALL_DOCDIR "${CMAKE_INSTALL_DATAROOTDIR}/doc/${PROJECT_NAME}")
endif()

mark_as_advanced(
  CMAKE_INSTALL_BINDIR
  CMAKE_INSTALL_LIBDIR
  CMAKE_INSTALL_INCLUDEDIR
  CMAKE_INSTALL_DATAROOTDIR
  CMAKE_INSTALL_DATADIR
  CMAKE_INSTALL_MANDIR
  CMAKE_INSTALL_DOCDIR
  )

set(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(FEELPP_PREFIX ${CMAKE_INSTALL_PREFIX})
if (NOT FEELPP_DATADIR )
  set(FEELPP_DATADIR ${CMAKE_INSTALL_PREFIX}/share/feelpp/feel )
endif()

#FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/applications/crb/templates/*")
#INSTALL(FILES ${files} DESTINATION share/feel/crb/templates COMPONENT Devel)

# documentation and examples
#  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/feel_get_tutorial.sh DESTINATION bin COMPONENT Doc)
#install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/manual/feelpp-manual.pdf DESTINATION share/doc/feel COMPONENT Doc)
IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/doc/api/html" )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/api/html DESTINATION share/doc/feelpp/feel COMPONENT Doc
    PATTERN ".svn" EXCLUDE PATTERN ".git" EXCLUDE )
ENDIF()


#FILE(GLOB files "quickstart/qs_*")
#INSTALL(FILES ${files} DESTINATION share/doc/feel/examples/ COMPONENT Doc)

#FILE(WRITE CMakeLists.txt.doc  "cmake_minimum_required(VERSION 2.8)
#Find_Package(Feel++ PATHS ${CMAKE_INSTALL_PREFIX}/share/feelpp/feel/cmake/modules/ )

#feelpp_add_application( qs_laplacian SRCS qs_laplacian.cpp INCLUDE_IN_ALL)
#feelpp_add_application( qs_stokes SRCS qs_stokes.cpp INCLUDE_IN_ALL)
#feelpp_add_application( qs_ns SRCS qs_ns.cpp INCLUDE_IN_ALL)
#")

#FILE(GLOB examples
#  "${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/tutorial/*.*pp")
#FILE(GLOB examplescfg
#  "${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/tutorial/*.cfg"
#  "${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/tutorial/*.geo" )

#INSTALL(FILES ${examples} DESTINATION share/doc/feel/examples/ COMPONENT Doc)
#INSTALL(FILES ${examplescfg} DESTINATION share/doc/feel/examples/ COMPONENT Doc)
#foreach(example ${examples} )
#  get_filename_component( EXAMPLE_TARGET_NAME ${example} NAME_WE )
#  get_filename_component( EXAMPLE_SRCS_NAME ${example} NAME )
#  FILE(APPEND CMakeLists.txt.doc "
# target feelpp_doc_${EXAMPLE_TARGET_NAME}
#feelpp_add_application( doc_${EXAMPLE_TARGET_NAME} SRCS ${EXAMPLE_SRCS_NAME} INCLUDE_IN_ALL)
#" )
#endforeach()
#foreach( example myapp mymesh myintegrals myfunctionspace mylaplacian mystokes)
#  FILE(APPEND CMakeLists.txt.doc "
#add_dependencies(tutorial feelpp_doc_${example})
#")
#endforeach()
#INSTALL(FILES CMakeLists.txt.doc DESTINATION share/doc/feel/examples/ COMPONENT Doc RENAME CMakeLists.txt)

#
# this target installs the libraries, header files and cmake files
#
set(_INSTALL_FEELPP_LIB_COMMAND ${CMAKE_COMMAND})

if ( FEELPP_HAS_GFLAGS )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/gflags/cmake_install.cmake")
endif()
if ( FEELPP_HAS_GLOG )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/glog/cmake_install.cmake")
endif()
if ( FEELPP_HAS_GINAC )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/ginac/cmake_install.cmake")
endif()
set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/eigen/cmake_install.cmake")
if(FEELPP_ENABLE_METIS)
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/metis/cmake_install.cmake")
endif()
if ( FEELPP_HAS_NLOPT )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/nlopt/cmake_install.cmake")
endif()
if ( FEELPP_HAS_IPOPT )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/ipopt/cmake_install.cmake")
endif()
if ( FEELPP_HAS_PYBIND11 ) #AND FEELPP_ENABLE_PYTHON AND FEELPP_ENABLE_PYTHON_WRAPPING )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/pybind11/cmake_install.cmake")
endif()

if ( FEELPP_HAS_MONGOCXX )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/mongocxx/src/cmake_install.cmake")
endif()

if ( TARGET feelpp_mesh_partitioner )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_BINARY_DIR}/applications/mesh/cmake_install.cmake")
endif()

if ( TARGET app-databases )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_BINARY_DIR}/applications/databases/cmake_install.cmake")
endif()

if ( TARGET feelpp_remotedata )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_BINARY_DIR}/tools/remotedata/cmake_install.cmake")
endif()

set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND}
  -DCMAKE_INSTALL_COMPONENT=Bin   -P "${CMAKE_BINARY_DIR}/tools/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/contrib/eigen/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Libs  -P "${CMAKE_BINARY_DIR}/feel/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/feel/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/cmake/modules/cmake_install.cmake"
  )

set( FEELPP_INSTALL_FEELPP_LIB_DEPENDS_TARGET contrib tools feelpp )
if ( TARGET feelpp_mesh_partitioner )
  list(APPEND FEELPP_INSTALL_FEELPP_LIB_DEPENDS_TARGET feelpp_mesh_partitioner )
endif()
if ( TARGET app-databases )
  list(APPEND FEELPP_INSTALL_FEELPP_LIB_DEPENDS_TARGET app-databases )
endif()
if ( TARGET feelpp_remotedata )
  list(APPEND FEELPP_INSTALL_FEELPP_LIB_DEPENDS_TARGET feelpp_remotedata )
endif()

add_custom_target(install-feelpp-lib
    DEPENDS ${FEELPP_INSTALL_FEELPP_LIB_DEPENDS_TARGET}
    COMMAND ${_INSTALL_FEELPP_LIB_COMMAND}
    )



if ( NOT TARGET install-feelpp-base )

  add_custom_target(install-feelpp-base
    DEPENDS
    install-feelpp-lib
    install-quickstart
    COMMAND
    "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=Geo
    -P "${CMAKE_BINARY_DIR}/cmake_install.cmake"
    "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=Quickstart
    -P "${CMAKE_BINARY_DIR}/quickstart/cmake_install.cmake"
    )
endif()

# install feel++ interpreter
if( 0 )#FEELPP_HAS_CLING_INTERPRETER )
    # We create the feel++ interpreter bash script in the binary dir.
    set( CLING_INSTALL_PREFIX ${CMAKE_BINARY_DIR} )
    include( ${CMAKE_SOURCE_DIR}/cmake/modules/feelpp.interpreter.bash.cmake )
    # We recreate the feel++ interpreter bash script for the cmake prefix.
    install( CODE
        "
        set( CLING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
        include( ${CMAKE_SOURCE_DIR}/cmake/modules/feelpp.interpreter.bash.cmake )
        "
        )
endif()
