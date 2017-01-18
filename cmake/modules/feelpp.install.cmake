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

if(FEELPP_ENABLE_METIS)
    set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND} -P "${CMAKE_BINARY_DIR}/contrib/metis/cmake_install.cmake")
endif()

set(_INSTALL_FEELPP_LIB_COMMAND ${_INSTALL_FEELPP_LIB_COMMAND}
  -P "${CMAKE_BINARY_DIR}/contrib/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Libs -P "${CMAKE_BINARY_DIR}/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Devel -P "${CMAKE_BINARY_DIR}/cmake_install.cmake")

if ( FEELPP_MINIMAL_BUILD )
  add_custom_target(install-feelpp-lib
    DEPENDS contrib feelpp 
    COMMAND ${_INSTALL_FEELPP_LIB_COMMAND}
    )
else()
  add_custom_target(install-feelpp-lib
    DEPENDS contrib feelpp feelpp_mesh_partitioner
    COMMAND ${_INSTALL_FEELPP_LIB_COMMAND}
    )
endif()

add_custom_target(install-libs-models-thermodyn
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  feelpp_model_thermodynamics
  COMMAND
      "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=Lib
      -P "${CMAKE_BINARY_DIR}/feel/feelmodels/thermodyn/cmake_install.cmake"
)

add_custom_target(install-apps-models-thermodyn
  DEPENDS
  install-libs-models-thermodyn
  feelpp_application_thermodyn_2d
  feelpp_application_thermodyn_3d
COMMAND
  "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=testcases
  -P "${CMAKE_BINARY_DIR}/applications/models/thermodyn/cmake_install.cmake"
  COMMAND
      "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=ModelApplications
      -P "${CMAKE_BINARY_DIR}/applications/models/thermodyn/cmake_install.cmake"
)
add_custom_target(install-libs-models-thermoelectric
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  feelpp_model_thermoelectric
  COMMAND
  "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=Lib
      -P "${CMAKE_BINARY_DIR}/feel/feelmodels/thermoelectric/cmake_install.cmake"
)

add_custom_target(install-apps-models-thermoelectric
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  install-libs-models-thermoelectric
  feelpp_application_thermoelectric_2d
  feelpp_application_thermoelectric_3d
    COMMAND
      "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=testcases
      -P "${CMAKE_BINARY_DIR}/applications/models/thermoelectric/cmake_install.cmake"
  COMMAND
      "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=ModelApplications
      -P "${CMAKE_BINARY_DIR}/applications/models/thermoelectric/cmake_install.cmake"
)

add_custom_target(install-apps-models-fluid
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  install-libs-models-thermodyn
  install-feelpp_model_fluidmec2dP2P1G1
  install-feelpp_model_fluidmec3dP2P1G1
  feelpp_application_fluid_2d
  feelpp_application_fluid_3d
  COMMAND
  "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=testcases
  -P "${CMAKE_BINARY_DIR}/applications/models/fluid/cmake_install.cmake"
  COMMAND
      "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=ModelApplications
      -P "${CMAKE_BINARY_DIR}/applications/models/fluid/cmake_install.cmake"
)

add_custom_target(install-apps-models-solid
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  install-feelpp_model_solidmec2dP1G1
  install-feelpp_model_solidmec2dP2G1
  install-feelpp_model_solidmec2dP1G1
  install-feelpp_model_solidmec3dP2G1
  feelpp_application_solid_2d
  feelpp_application_solid_3d
  feelpp_application_stress_2d
  feelpp_application_stress_3d
  feelpp_application_solenoid_3d
  COMMAND
  "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=testcases
  -P "${CMAKE_BINARY_DIR}/applications/models/solid/cmake_install.cmake"
  COMMAND
      "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=ModelApplications
      -P "${CMAKE_BINARY_DIR}/applications/models/solid/cmake_install.cmake"
)

add_custom_target(install-apps-models-fsi
  DEPENDS
  install-feelpp-lib
  install-feelpp-models-common
  install-feelpp_model_fsi_2dP2P1G1_2dP1G1
  install-feelpp_model_fsi_3dP2P1G1_3dP1G1
  install-feelpp_model_fsi_3dP2P1G2_3dP2G2
  feelpp_application_fsi_2d
  feelpp_application_fsi_3d
  feelpp_application_fsi_3d_g2
  COMMAND
  "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=testcases
  -P "${CMAKE_BINARY_DIR}/applications/models/fsi/cmake_install.cmake"
  COMMAND
      "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=ModelApplications
      -P "${CMAKE_BINARY_DIR}/applications/models/fsi/cmake_install.cmake"
)

if ( NOT TARGET install-feelpp-base )
  #
  # this target installs the libraries, header files, cmake files and sample applications
  #
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
