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
set(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(FEELPP_PREFIX ${CMAKE_INSTALL_PREFIX})
if (NOT FEELPP_DATADIR )
  set(FEELPP_DATADIR ${CMAKE_INSTALL_PREFIX}/share/feel )
endif()


FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/applications/crb/templates/*")
INSTALL(FILES ${files} DESTINATION share/feel/crb/templates COMPONENT Devel)

# documentation and examples
#  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/feel_get_tutorial.sh DESTINATION bin COMPONENT Doc)
#install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/manual/feelpp-manual.pdf DESTINATION share/doc/feel COMPONENT Doc)
IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/doc/api/html" )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/api/html DESTINATION share/doc/feel COMPONENT Doc
    PATTERN ".svn" EXCLUDE PATTERN ".git" EXCLUDE )
ENDIF()

FILE(GLOB files "quickstart/qs_*")
INSTALL(FILES ${files} DESTINATION share/doc/feel/examples/ COMPONENT Doc)

FILE(WRITE CMakeLists.txt.doc  "cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH \"${CMAKE_INSTALL_PREFIX}/share/feel/cmake/modules/\")
Find_Package(Feel++)


feelpp_add_application( qs_laplacian SRCS qs_laplacian.cpp INCLUDE_IN_ALL)
feelpp_add_application( qs_stokes SRCS qs_stokes.cpp INCLUDE_IN_ALL)
feelpp_add_application( qs_ns SRCS qs_ns.cpp INCLUDE_IN_ALL)
")

#FILE(GLOB examples
#  "${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/tutorial/*.*pp")
#FILE(GLOB examplescfg
#  "${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/tutorial/*.cfg"
#  "${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/tutorial/*.geo" )

INSTALL(FILES ${examples} DESTINATION share/doc/feel/examples/ COMPONENT Doc)
INSTALL(FILES ${examplescfg} DESTINATION share/doc/feel/examples/ COMPONENT Doc)
foreach(example ${examples} )
  get_filename_component( EXAMPLE_TARGET_NAME ${example} NAME_WE )
  get_filename_component( EXAMPLE_SRCS_NAME ${example} NAME )
  FILE(APPEND CMakeLists.txt.doc "
# target feelpp_doc_${EXAMPLE_TARGET_NAME}
feelpp_add_application( doc_${EXAMPLE_TARGET_NAME} SRCS ${EXAMPLE_SRCS_NAME} INCLUDE_IN_ALL)
" )
endforeach()
#foreach( example myapp mymesh myintegrals myfunctionspace mylaplacian mystokes)
#  FILE(APPEND CMakeLists.txt.doc "
#add_dependencies(tutorial feelpp_doc_${example})
#")
#endforeach()
INSTALL(FILES CMakeLists.txt.doc DESTINATION share/doc/feel/examples/ COMPONENT Doc RENAME CMakeLists.txt)

#
# this target installs the libraries, header files and cmake files
#
set(_INSTALL_FEELPP_LIB_COMMAND ${CMAKE_COMMAND})

set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND}
  -P "${CMAKE_SOURCE_DIR}/cmake/modules/feelpp.install.config.cmake")

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

if ( TARGET feelpp_mesh_partitioner )
  set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_BINARY_DIR}/applications/mesh/cmake_install.cmake")
endif()

set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND}
  -P "${CMAKE_BINARY_DIR}/tools/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/feel/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/feel/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/cmake/modules/cmake_install.cmake"
  )

set( FEELPP_INSTALL_FEELPP_LIB_DEPENDS_TARGET contrib tools feelpp )
if ( TARGET feelpp_mesh_partitioner )
  list(APPEND FEELPP_INSTALL_FEELPP_LIB_DEPENDS_TARGET feelpp_mesh_partitioner )
endif()
  add_custom_target(install-feelpp-lib
    DEPENDS ${FEELPP_INSTALL_FEELPP_LIB_DEPENDS_TARGET}
    COMMAND ${_INSTALL_FEELPP_LIB_COMMAND}
    )


add_custom_target(install-feelpp-models-common
  DEPENDS
  feelpp-models-common
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelcore/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelcore/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelalg/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelalg/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelmesh/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelmesh/cmake_install.cmake"
)

add_custom_target(install-libs-models-thermodyn
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  feelpp_model_thermodynamics
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/thermodyn/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/thermodyn/cmake_install.cmake"
)

add_custom_target(install-apps-models-thermodyn
  DEPENDS
  install-libs-models-thermodyn
  feelpp_toolbox_thermodyn_2d
  feelpp_toolbox_thermodyn_3d
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=testcases -P "${CMAKE_BINARY_DIR}/toolboxes/thermodyn/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_BINARY_DIR}/toolboxes/thermodyn/cmake_install.cmake"
)
add_custom_target(install-libs-models-thermoelectric
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  feelpp_model_thermoelectric
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/thermoelectric/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/thermoelectric/cmake_install.cmake"
)

add_custom_target(install-apps-models-thermoelectric
  DEPENDS
  install-libs-models-thermoelectric
  feelpp_toolbox_thermoelectric_2d
  feelpp_toolbox_thermoelectric_3d
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=testcases -P "${CMAKE_BINARY_DIR}/toolboxes/thermoelectric/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_BINARY_DIR}/toolboxes/thermoelectric/cmake_install.cmake"
)

add_custom_target(install-libs-models-fluid
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  install-libs-models-thermodyn
  feelpp_model_fluidmechanics
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fluid/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fluid/cmake_install.cmake"
)
add_custom_target(install-apps-models-fluid
  DEPENDS
  install-libs-models-fluid
  feelpp_toolbox_fluid_2d
  feelpp_toolbox_fluid_3d
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=testcases -P "${CMAKE_BINARY_DIR}/toolboxes/fluid/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_BINARY_DIR}/toolboxes/fluid/cmake_install.cmake"
)

add_custom_target(install-libs-models-solid
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  feelpp_model_solidmechanics
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/solid/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/solid/cmake_install.cmake"
)
add_custom_target(install-apps-models-solid
  DEPENDS
  install-libs-models-solid
  feelpp_toolbox_solid_2d
  feelpp_toolbox_solid_3d
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=testcases -P "${CMAKE_BINARY_DIR}/toolboxes/solid/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_BINARY_DIR}/toolboxes/solid/cmake_install.cmake"
)

add_custom_target(install-libs-models-fsi
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  feelpp_model_fsi
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fsi/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fsi/cmake_install.cmake"
)
add_custom_target(install-apps-models-fsi
  DEPENDS
  install-libs-models-fsi
  feelpp_toolbox_fsi_2d
  feelpp_toolbox_fsi_3d
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=testcases -P "${CMAKE_BINARY_DIR}/toolboxes/fsi/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_BINARY_DIR}/toolboxes/fsi/cmake_install.cmake"
)

if ( NOT TARGET install-feelpp-base )
  #
  # this target installs the libraries, header files, cmake files and sample applications
  #
  if ( FEELPP_HAS_ACUSIM )
    set(FEELPP_DATABASES_APPS  feelpp_databases_converter_acusim     feelpp_databases_export     feelpp_databases_pod)
  else()
    set(FEELPP_DATABASES_APPS  feelpp_databases_export     feelpp_databases_pod)
  endif()
  MESSAGE(STATUS "toto ${FEELPP_DATABASES_APPS}")

  add_custom_target(install-feelpp-base
    DEPENDS
    install-feelpp-lib
    install-quickstart
    install-app-databases
    COMMAND
    "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=Geo
    -P "${CMAKE_BINARY_DIR}/cmake_install.cmake"
    "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=Quickstart
    -P "${CMAKE_BINARY_DIR}/quickstart/cmake_install.cmake"
    "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=Databases
    -P "${CMAKE_BINARY_DIR}/applications/databases/cmake_install.cmake"
    )
endif()

#
# Generic install target for feel++
#
add_custom_target(install-feelpp
  DEPENDS
  install-feelpp-base
)

# install feel++ applications
add_custom_target(install-feelpp-apps
  DEPENDS
  install-feelpp-base
  install-apps-models-thermoelectric
  install-apps-models-thermodyn
  install-apps-models-fluid
  install-apps-models-solid
  install-apps-models-fsi
)

# install feel++ interpreter
if( FEELPP_ENABLE_INTERPRETER )
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
