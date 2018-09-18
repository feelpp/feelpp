###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 01 Sep 2018
#
#  Copyright (C) 2018 Feel++ Consortium
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
add_library(feelpp_contrib INTERFACE)

set( FEELPP_HAS_GFLAGS 1 )
add_subdirectory(gflags)
if ( FEELPP_HAS_GFLAGS )
  set(GFLAGS_IS_SUBPROJECT TRUE)
  set(GFLAGS_NAMESPACE "google;gflags")
  #add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gflags)
  list(INSERT FEELPP_LIBRARIES 0 feelpp_gflags)

  #add_dependencies(contrib feelpp_gflags feelpp_gflags_shared feelpp_gflags_nothreads_shared)
  #target_include_directories(feelpp_gflags_nothreads_shared BEFORE PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/gflags/include ${CMAKE_CURRENT_BINARY_DIR}/gflags/include/gflags)

  target_include_directories(feelpp_contrib
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/feelpp/contrib/gflags/include>
    $<INSTALL_INTERFACE:include/feelpp>
    )
  target_link_libraries(feelpp_contrib INTERFACE feelpp_gflags)

  set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} GFlags/Contrib")
  
  # for GLog
  set(gflags_LIBRARIES feelpp_gflags)
  set(gflags_FOUND 1)
endif()

set( FEELPP_HAS_GLOG 1 )
add_subdirectory(glog)
if ( FEELPP_HAS_GLOG )
  target_include_directories(feelpp_contrib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/feelpp/contrib/glog/>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/feelpp/contrib/glog/src>
    $<INSTALL_INTERFACE:include/feelpp>
    )
  target_link_libraries(feelpp_contrib INTERFACE feelpp_glog)
    #     INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/glog/ ${CMAKE_CURRENT_SOURCE_DIR}/glog/src)
  #add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/glog)
  list(INSERT FEELPP_LIBRARIES 0 feelpp_glog)
  #add_dependencies(contrib feelpp_glog)
  set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} GLog/Contrib")
endif()

set( FEELPP_HAS_GINAC 1 )
add_subdirectory(ginac)
if ( FEELPP_HAS_GINAC )
  #
  # cln and ginac
  #
  find_package(CLN)

  target_compile_definitions(feelpp_contrib INTERFACE -DIN_GINAC -DHAVE_LIBDL)

  #link_directories( ${CMAKE_BINARY_DIR}/contrib/ginac/ginac)

  target_include_directories( feelpp_contrib BEFORE
    INTERFACE
    ${CLN_INCLUDE_DIR}
    $<BUILD_INTERFACE:${FEELPP_SOURCE_DIR}/contrib/ginac/>
    $<BUILD_INTERFACE:${FEELPP_BUILD_DIR}/ginac/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/contrib/ginac/ginac>
    $<BUILD_INTERFACE:${FEELPP_BUILD_DIR}/contrib/ginac/ginac>
    $<INSTALL_INTERFACE:include/feelpp>
    $<INSTALL_INTERFACE:include/feelpp/ginac>
    )
  set(DL_LIBS ${CMAKE_DL_LIBS})

  target_link_libraries(feelpp_contrib INTERFACE feelpp_ginac)
  #add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ginac)
  #add_dependencies(contrib feelpp_ginac ginsh)
  list(INSERT FEELPP_LIBRARIES 0 feelpp_ginac)
endif()


#
# Eigen
#
add_subdirectory(eigen)
if ( FEELPP_ENABLE_SYSTEM_EIGEN3 )
  FIND_PACKAGE(Eigen3)
  MESSAGE(STATUS "Eigen3 system found:")
  MESSAGE("EIGEN_INCLUDE_DIR=${EIGEN_INCLUDE_DIR}")
  MESSAGE("EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}")
  MESSAGE(STATUS "Adding unsupported headers to EIGEN3_INCLUDE_DIR:")
  set( EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/unsupported)
  MESSAGE("EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}")
endif()
if (NOT EIGEN3_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/feelpp/feel AND EXISTS ${CMAKE_SOURCE_DIR}/feelpp/contrib )
  option(EIGEN_BUILD_PKGCONFIG "Build pkg-config .pc file for Eigen" OFF)
  
  set( EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/feelpp/contrib/eigen ${CMAKE_SOURCE_DIR}/feelpp/contrib/eigen/unsupported )

  SET(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Eigen3/Contrib" )
elseif( EIGEN3_FOUND )
  SET(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Eigen3/System" )
else()
  find_path(EIGEN3_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
    PATHS
    $ENV{FEELPP_DIR}/include/feelpp
    NO_DEFAULT_PATH
    )
endif()

target_link_libraries(feelpp_contrib INTERFACE eigen )
get_property(INCD TARGET eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES )
message(STATUS "incd1: ${INCD}")
get_property(INCD TARGET feelpp_contrib PROPERTY INTERFACE_INCLUDE_DIRECTORIES )
message(STATUS "incd2: ${INCD}")

SET(FEELPP_HAS_EIGEN3 1)
if ( FEELPP_HAS_EIGEN3 )
  #add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/contrib/eigen)
  unset(INCLUDE_INSTALL_DIR CACHE)
  unset(CMAKEPACKAGE_INSTALL_DIR CACHE)
  unset(PKGCONFIG_INSTALL_DIR CACHE)
endif()
message(STATUS "[feelpp] eigen3 headers: ${EIGEN3_INCLUDE_DIR}" )


#FIND_PACKAGE(Eigen2 REQUIRED)
#INCLUDE_DIRECTORIES( ${Eigen2_INCLUDE_DIR} )
#add_subdirectory(contrib/eigen)
#INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/contrib/eigen )
#add_definitions( -DEIGEN_NO_STATIC_ASSERT )

#
# METIS
#
option( FEELPP_ENABLE_METIS "Enable Metis Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if(FEELPP_ENABLE_METIS)
  add_subdirectory(metis)
  SET(FEELPP_HAS_METIS 1)
  TARGET_INCLUDE_DIRECTORIES(feelpp_contrib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/feelpp/contrib/metis/include>
    $<INSTALL_INTERFACE:include/feelpp>)
  target_link_libraries(feelpp_contrib INTERFACE feelpp_metis)
  SET(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Metis/Contrib" )
endif( FEELPP_ENABLE_METIS )


if ( NOT FEELPP_HAS_PARMETIS )
  option( FEELPP_ENABLE_PARMETIS "Enable Parmetis Support" OFF )
  if(FEELPP_ENABLE_PARMETIS)
    FIND_LIBRARY(PARMETIS_LIBRARY
      NAMES
      parmetis
      PATHS
      $ENV{PETSC_DIR}/lib
      $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
      )
    IF( PARMETIS_LIBRARY )
      message(STATUS "[feelpp] Parmetis: ${PARMETIS_LIBRARY}" )
      SET(FEELPP_LIBRARIES ${PARMETIS_LIBRARY} ${FEELPP_LIBRARIES})
      set(FEELPP_HAS_PARMETIS 1)
    ENDIF()
  endif()
endif()




option( FEELPP_ENABLE_NLOPT "Enable NLOPT (NonLinear Optimisation Library)" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )

if ( FEELPP_ENABLE_NLOPT )
  feelppContribPrepare( nlopt )

  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    add_subdirectory( nlopt )
    set(NLOPT_INCLUDE_DIR ${FEELPP_SOURCE_DIR}/contrib/nlopt/api  ${FEELPP_BINARY_DIR}/contrib/nlopt/api)
    
    target_link_libraries(feelpp_contrib INTERFACE feelpp_nlopt )
    SET( FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} NLOpt/Contrib")
    SET(FEELPP_HAS_NLOPT 1)
    list(INSERT FEELPP_LIBRARIES 0 feelpp_nlopt)
  endif( FEELPP_CONTRIB_PREPARE_SUCCEED )
endif( FEELPP_ENABLE_NLOPT )

if( FEELPP_HAS_IPOPT )
  add_subdirectory(ipopt)
endif()
## No cmake support for HPDDM!
#if( FEELPP_HAS_HPDDM )
#  add_subdirectory(hpddm)
#endif()
option( FEELPP_ENABLE_PYBIND11 "Enable PYBIND11" ON )
if ( FEELPP_ENABLE_PYBIND11 )
  feelppContribPrepare( pybind11 )
  if ( FEELPP_CONTRIB_PREPARE_SUCCEED )
    add_subdirectory(pybind11)
    set( FEELPP_HAS_PYBIND11 1 )
    target_link_libraries(feelpp_contrib INTERFACE pybind11 )
    target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_PYBIND11=1 )
    target_include_directories(feelpp_contrib INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/pybind11/include>
      $<INSTALL_INTERFACE:include/feelpp/pybind11> )
    set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} PyBind11/Contrib")
  endif()
endif()

#
# MONGOCXX
#
option( FEELPP_ENABLE_MONGOCXX "Enable Mongocxx" OFF )
if ( FEELPP_ENABLE_MONGOCXX )
  feelppContribPrepare( mongocxx )
  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    add_subdirectory(mongocxx)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${FEELPP_SOURCE_DIR}/contrib/mongocxx/cmake")
    find_package(LibBSON 1.5 )
    find_package(LibMongoC 1.5 )
    if ( LIBBSON_FOUND AND LIBMONGOC_FOUND )
      message(STATUS "[feelpp] found LibBSON and LibMongoC")
      set(BSONCXX_INLINE_NAMESPACE "")
      set(BSONCXX_HEADER_INSTALL_DIR "include/feelpp/bsoncxx/" CACHE INTERNAL "")

      include_directories(${LIBBSON_INCLUDE_DIRS} ${LIBMONGOC_INCLUDE_DIRS})
      set(MONGOCXX_LIBRARIES "mongocxx_static ${LIBMONGOC_LIBRARIES} ${LIBBSON_LIBRARIES}")
      set(MONGOCXX_INCLUDE_DIRS "${LIBBSON_INCLUDE_DIRS} ${LIBMONGOC_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/src     ${CMAKE_CURRENT_BINARY_DIR}/src")
      set(LIBBSON_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
      list(APPEND LIBBSONCXX_INCLUDE_DIRS ${FEELPP_SOURCE_DIR}/contrib/mongocxx/src ${FEELPP_BINARY_DIR}/contrib/mongocxx/src)
      list(APPEND LIBMONGOCXX_INCLUDE_DIRS ${FEELPP_SOURCE_DIR}/contrib/mongocxx/src ${FEELPP_BINARY_DIR}/contrib/mongocxx/src)
      include_directories(${LIBMONGOCXX_INCLUDE_DIRS} ${LIBBSONCXX_INCLUDE_DIRS})
      set(FEELPP_HAS_BSONCXX 1)
      set(FEELPP_HAS_MONGOCXX 1)
      set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} MongoCXX/Contrib" )
      list(APPEND MONGO_LIBRARIES feelpp_mongocxx feelpp_bsoncxx ${LIBMONGOC_LIBRARIES} ${LIBBSON_LIBRARIES} )
      list(INSERT FEELPP_LIBRARIES 0 ${MONGO_LIBRARIES})
    endif()
  else()
    message(WARNING "MongoCXX was not found on your system. Either install it or set FEELPP_ENABLE_MONGOCXX to OFF.")
  endif()
endif( FEELPP_ENABLE_MONGOCXX )


add_subdirectory(gmsh)

#include(feelpp.module.mongocxx)
#include(feelpp.module.hpddm)
#include(feelpp.module.nlopt)
#include(feelpp.module.ipopt)
#include(feelpp.module.cereal)
#include(feelpp.module.paralution)
#include(feelpp.module.jsonlab)

# Add an info message to be displayed at the end of the cmake process.
if( FEELPP_CONTRIB_SUBMODULE_UPDATED )
  list( APPEND FEELPP_MESSAGE_INFO_END "Feel++ submodules already initialized!\nPlease make sure submodules are up to date (run `git submodule update --init --recursive` in source directory)" )
  set( FEELPP_MESSAGE_INFO_END ${FEELPP_MESSAGE_INFO_END} )
endif()


INSTALL(TARGETS feelpp_contrib DESTINATION lib/ COMPONENT Libs EXPORT feelpp-contrib-export-targets)
install(EXPORT feelpp-contrib-export-targets DESTINATION share/feelpp/feel/cmake/modules)

set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB}" PARENT_SCOPE)
