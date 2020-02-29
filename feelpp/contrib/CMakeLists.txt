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
add_library( Feelpp::feelpp_contrib ALIAS feelpp_contrib)
set( FEELPP_HAS_GFLAGS 1 )
set( FEELPP_HAS_GFLAGS 1 PARENT_SCOPE )
set( GFLAGS_IS_SUBPROJECT TRUE )
add_subdirectory(gflags)
if ( FEELPP_HAS_GFLAGS )
  set(GFLAGS_IS_SUBPROJECT TRUE)
  set(GFLAGS_NAMESPACE "google;gflags")
  #add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gflags)
  list(INSERT FEELPP_LIBRARIES 0 feelpp_gflags)

  #add_dependencies(contrib feelpp_gflags feelpp_gflags_shared feelpp_gflags_nothreads_shared)
  #target_include_directories(feelpp_gflags_nothreads_shared BEFORE PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/gflags/include ${CMAKE_CURRENT_BINARY_DIR}/gflags/include/gflags)

  target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_GFLAGS )
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
set( FEELPP_HAS_GLOG 1 PARENT_SCOPE )
add_subdirectory(glog)
if ( FEELPP_HAS_GLOG )
  target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_GLOG )
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

add_subdirectory(ginac)
#
# cln and ginac
#
find_package(CLN)

target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_GINAC IN_GINAC HAVE_LIBDL)

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



#
# Eigen
#
if ( FEELPP_ENABLE_SYSTEM_EIGEN3 )
  FIND_PACKAGE(Eigen3 REQUIRED)
  MESSAGE(STATUS "Eigen3 system found:")
else()

  set(EIGEN_LEAVE_TEST_IN_ALL_TARGET 1)
  set( EIGEN_INCLUDE_INSTALL_DIR "include/feelpp/eigen" )
  add_subdirectory(eigen)

  if (NOT EIGEN3_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/feelpp/feel AND EXISTS ${CMAKE_SOURCE_DIR}/feelpp/contrib )
    option(EIGEN_BUILD_PKGCONFIG "Build pkg-config .pc file for Eigen" OFF)
    
    set( EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/feelpp/contrib/eigen ${CMAKE_SOURCE_DIR}/feelpp/contrib/eigen/unsupported )
    target_include_directories( feelpp_contrib BEFORE
      INTERFACE
      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/contrib/eigen/>
      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/feelpp/contrib/eigen/unsupported>
      $<INSTALL_INTERFACE:include/feelpp/eigen>
      $<INSTALL_INTERFACE:include/feelpp/eigen/unsupported>
      )

    SET(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Eigen3/Contrib" )
  else()
    message(FATAL_ERROR "Eigen contrib not found" )
    # find_path(EIGEN3_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
    #   PATHS
    #   $ENV{FEELPP_DIR}/include/feelpp
    #   NO_DEFAULT_PATH
    #   )
  endif()

endif()

# there is an issue with eigen+mkl
if ( 0 ) #MKL_FOUND )
  target_compile_definitions(eigen INTERFACE EIGEN_USE_MKL_ALL)
  target_link_libraries(eigen INTERFACE ${MKL_LIBRARIES})  # error with cmake < 3.13, see CMP0079, use instead feelpp_contrib
  #target_link_libraries(feelpp_contrib INTERFACE ${MKL_LIBRARIES})
  target_include_directories(eigen INTERFACE ${MKL_INCLUDE_DIRS})
  target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_MKL )
endif()
target_compile_definitions(feelpp_contrib INTERFACE FEELP_HAS_EIGEN3 )
target_link_libraries(feelpp_contrib INTERFACE eigen )

SET(FEELPP_HAS_EIGEN3 1)
SET(FEELPP_HAS_EIGEN3 1 PARENT_SCOPE)
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
  #SET(FEELPP_HAS_METIS 1)
  #SET(FEELPP_HAS_METIS 1 PARENT_SCOPE)
  target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_METIS )
  TARGET_INCLUDE_DIRECTORIES(feelpp_contrib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/feelpp/contrib/metis/include>
    $<INSTALL_INTERFACE:include/feelpp/metis>)
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
      set(FEELPP_HAS_PARMETIS 1 PARENT_SCOPE)
    ENDIF()
  endif()
endif()




option( FEELPP_ENABLE_NLOPT "Enable NLOPT (NonLinear Optimisation Library)" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )

if ( FEELPP_ENABLE_NLOPT )
  feelppContribPrepare( nlopt )

  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    add_subdirectory( nlopt )
    set(NLOPT_INCLUDE_DIR ${FEELPP_SOURCE_DIR}/contrib/nlopt/api  ${FEELPP_BINARY_DIR}/contrib/nlopt/api)
    target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_NLOPT )
    target_link_libraries(feelpp_contrib INTERFACE feelpp_nlopt )
    SET( FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} NLOpt/Contrib")
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
set(PYBIND11_CPP_STANDARD -std=c++${FEELPP_STD_CPP})
if ( FEELPP_ENABLE_PYBIND11 )
  feelppContribPrepare( pybind11 )
  if ( FEELPP_CONTRIB_PREPARE_SUCCEED )
    add_subdirectory(pybind11)
    target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_PYBIND11 )
    target_link_libraries(feelpp_contrib INTERFACE pybind11 )
    #target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_PYBIND11=1 )
    target_include_directories(feelpp_contrib INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/pybind11/include>
      $<INSTALL_INTERFACE:include/feelpp/pybind11> )
    set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} PyBind11/Contrib")
  endif()
endif()


# MongoCXX
option( FEELPP_ENABLE_MONGOCXX "Enable Mongocxx" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if ( FEELPP_ENABLE_MONGOCXX )
  feelppContribPrepare( mongocxx )
  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/mongocxx/cmake)
    set(LIBBSON_REQUIRED_VERSION 1.5.0)
    set(LIBMONGOC_REQUIRED_VERSION 1.5.0)
    find_package(LibBSON ${LIBBSON_REQUIRED_VERSION})
    find_package(LibMongoC ${LIBMONGOC_REQUIRED_VERSION} )
    if ( LIBBSON_FOUND AND LIBMONGOC_FOUND )
      message(STATUS "[feelpp] found LibBSON and LibMongoC")
      add_subdirectory(mongocxx)
      target_link_libraries(feelpp_contrib INTERFACE feelpp_mongocxx)

      if ( 0 )
      set(MONGOCXX_LIBRARIES "mongocxx_static ${LIBMONGOC_LIBRARIES} ${LIBBSON_LIBRARIES}")
      set(MONGOCXX_INCLUDE_DIRS "${LIBBSON_INCLUDE_DIRS} ${LIBMONGOC_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/src     ${CMAKE_CURRENT_BINARY_DIR}/src")
      endif()
      target_include_directories(feelpp_contrib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mongocxx/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mongocxx/src>
        $<INSTALL_INTERFACE:include/feelpp>
        )
      target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_BSONCXX FEELPP_HAS_MONGOCXX )
    else()
        if( NOT LIBBSON_FOUND )
            message(STATUS "[feelpp] BSON library was not found on your system!")
        endif()
        if( NOT LIBBSON_FOUND )
            message(STATUS "[feelpp] MONGODB library was not found on your system!")
        endif()
        message(STATUS "[feelpp] MONGOCXX dependencies missing. Either install it or set FEELPP_ENABLE_MONGOCXX to OFF for suppress this message.")
    endif()
  endif()
endif()


option( FEELPP_ENABLE_KWSYS "Enable KWSys" ON )
# Kwsys configure
set(KWSYS_NAMESPACE feelpp_kwsys CACHE INTERNAL "")
#  # Enable all components.
#  set(KWSYS_ENABLE_C 1)
#  set(KWSYS_USE_Base64 1)
#  set(KWSYS_USE_Directory 1)
#  set(KWSYS_USE_DynamicLoader 1)
#  set(KWSYS_USE_Encoding 1)
#  set(KWSYS_USE_Glob 1)
#  set(KWSYS_USE_MD5 1)
#  set(KWSYS_USE_Process 1)
#  set(KWSYS_USE_RegularExpression 1)
set(KWSYS_USE_System 1 CACHE INTERNAL "")
#  set(KWSYS_USE_SystemTools 1)
#  set(KWSYS_USE_CommandLineArguments 1)
#  set(KWSYS_USE_Terminal 1)
#  set(KWSYS_USE_IOStream 1)
#  set(KWSYS_USE_FStream 1)
#  set(KWSYS_USE_String 1)
set(KWSYS_USE_SystemInformation 1 CACHE INTERNAL "")
#  set(KWSYS_USE_ConsoleBuf 1)

# install rules.
set(KWSYS_INSTALL_BIN_DIR ${CMAKE_INSTALL_BINDIR} CACHE INTERNAL "")
set(KWSYS_INSTALL_LIB_DIR ${CMAKE_INSTALL_LIBDIR} CACHE INTERNAL "")
set(KWSYS_INSTALL_INCLUDE_DIR include/feelpp CACHE INTERNAL "")
set(KWSYS_INSTALL_COMPONENT_NAME_RUNTIME Runtime CACHE INTERNAL "")
set(KWSYS_INSTALL_COMPONENT_NAME_DEVELOPMENT Development CACHE INTERNAL "")
set(KWSYS_INSTALL_EXPORT_NAME feelpp-contrib-export-targets CACHE INTERNAL "")
if ( FEELPP_ENABLE_KWSYS )
  feelppContribPrepare( kwsys )
  if ( FEELPP_CONTRIB_PREPARE_SUCCEED )
    add_subdirectory(kwsys)
    target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_KWSYS )
    target_link_libraries(feelpp_contrib INTERFACE feelpp_kwsys )
    target_include_directories(feelpp_contrib INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/kwsys/include>
    $<INSTALL_INTERFACE:include/feelpp/feelpp_kwsys> )
    set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} KWSys/Contrib")
  endif()
endif()


option( FEELPP_ENABLE_FMILIB "Enable FMLIB C API" ON )
# FMILib configure
# install rules.
set(FMILIB_BUILD_STATIC_LIB OFF)
set(FMILIB_BUILD_SHARED_LIB ON)
set(FMILIB_BUILD_TESTS OFF)
set(FMILIB_INSTALL_SUBLIBS ON)
set(FMILIB_INSTALL_BIN_DIR ${CMAKE_INSTALL_BINDIR} CACHE INTERNAL "")
set(FMILIB_INSTALL_LIB_DIR ${CMAKE_INSTALL_LIBDIR} CACHE INTERNAL "")
set(FMILIB_INSTALL_INCLUDE_DIR include/feelpp/fmilib CACHE INTERNAL "")
set(FMILIB_INSTALL_COMPONENT_NAME_RUNTIME Runtime CACHE INTERNAL "")
set(FMILIB_INSTALL_COMPONENT_NAME_DEVELOPMENT Development CACHE INTERNAL "")
set(FMILIB_INSTALL_EXPORT_NAME feelpp-contrib-export-targets CACHE INTERNAL "")
if ( FEELPP_ENABLE_FMILIB )
  feelppContribPrepare( fmilib )
  if ( FEELPP_CONTRIB_PREPARE_SUCCEED )
    add_subdirectory(fmilib)
    target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_FMILIB )
    target_link_libraries(feelpp_contrib INTERFACE feelpp_fmilib )
    target_include_directories(feelpp_contrib INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/fmilib/>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/contrib/fmilib/>
      $<INSTALL_INTERFACE:include/feelpp/fmilib> )
    set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} FmiLib/Contrib")
  endif()
endif()

if ( FEELPP_HAS_GMSH )
  add_subdirectory(gmsh)
  if( FEELPP_HAS_GMSH_LIBRARY AND ( NOT FEELPP_HAS_GMSH_API OR FEELPP_GMSH_API_MISSING_HEADERS ) )
    if ( TARGET feelpp_gmsh_get_files )
      add_dependencies( feelpp_contrib feelpp_gmsh_get_files )
    endif()
    if ( DEFINED FEELPP_GMSH_MISSING_HEADERS_DIR )
      target_include_directories(feelpp_contrib INTERFACE
        $<BUILD_INTERFACE:${FEELPP_GMSH_MISSING_HEADERS_DIR}>
        $<INSTALL_INTERFACE:include/feelpp/gmsh> )
    endif()
  endif()
endif()

feelppContribPrepare( tabulate )
if( FEELPP_CONTRIB_PREPARE_SUCCEED )
  add_subdirectory( tabulate )
  target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_TABULATE )
  target_link_libraries(feelpp_contrib INTERFACE tabulate )
  target_include_directories(feelpp_contrib INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/tabulate/include>
    $<INSTALL_INTERFACE:include/feelpp> )
  set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Tabulate/Contrib")
endif()
feelppContribPrepare( indicators )
if( FEELPP_CONTRIB_PREPARE_SUCCEED )
  add_subdirectory( indicators )
  target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_INDICATORS )
  target_link_libraries(feelpp_contrib INTERFACE indicators )
  target_include_directories(feelpp_contrib INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/indicators/include>
    $<INSTALL_INTERFACE:include/feelpp> )
  set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Indicators/Contrib")
endif()
  
add_subdirectory( boost )
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
install(EXPORT feelpp-contrib-export-targets NAMESPACE Feelpp:: DESTINATION share/feelpp/feel/cmake/modules)

