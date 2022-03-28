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
message(STATUS "============================================================")
message(STATUS "Contrib libraries")
include(FetchContent)
add_library(feelpp_contrib INTERFACE)
add_library( Feelpp::feelpp_contrib ALIAS feelpp_contrib)

set( CMAKE_INSTALL_DATADIR "share/feelpp")
if ( NOT ${CMAKE_INSTALL_INCLUDEDIR} MATCHES "feelpp$")
    set( CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/feelpp")
endif()

message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring NAPP")
set(FETCHCONTENT_SOURCE_DIR_NAPP ${CMAKE_CURRENT_SOURCE_DIR}/napp)
FetchContent_Declare(napp GIT_REPOSITORY https://github.com/vincentchabannes/napp.git GIT_TAG 1cba0026b9ea7581902017074071a0da1a2e7a94) #v0.1.0 )
set(NAPP_CONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/feel/cmake/modules")
FetchContent_MakeAvailable(napp)
target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_NAPP )
target_link_libraries( feelpp_contrib INTERFACE napp )
install(TARGETS napp EXPORT feelpp-contrib-export-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring FMT")
set(FETCHCONTENT_SOURCE_DIR_FMT ${CMAKE_CURRENT_SOURCE_DIR}/fmt)
FetchContent_Declare(fmt GIT_REPOSITORY https://github.com/feelpp/fmt.git GIT_TAG b6f4ceaed0a0a24ccf575fab6c56dd50ccf6f1a9) 

set(FMT_INSTALL ON)
set(FMT_CMAKE_DIR "${CMAKE_INSTALL_DATADIR}/feel/cmake/modules")
FetchContent_MakeAvailable(fmt)
set_target_properties(fmt PROPERTIES OUTPUT_NAME "feelpp_fmt")
target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_FMT )
target_include_directories(fmt INTERFACE
    $<BUILD_INTERFACE:${fmt_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include/feelpp> )
target_link_libraries( feelpp_contrib INTERFACE fmt )
install(TARGETS fmt EXPORT feelpp-contrib-export-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring GFlags")
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

message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring GLog")
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

message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring Ginac")
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


message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring MKL")
option( FEELPP_ENABLE_MKL "Enable MKL Support" ON )
if ( FEELPP_ENABLE_MKL )
  find_package(MKL)
  target_include_directories( feelpp_contrib BEFORE INTERFACE ${MKL_INCLUDE_DIRS})
  target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_MKL )
endif()


# Eigen
#


include(FetchContent)
set(EIGEN_BUILD_DOC OFF)
set(BUILD_TESTING OFF)
set(EIGEN_LEAVE_TEST_IN_ALL_TARGET 0)
set(EIGEN_TEST_NOQT 1)
set(FETCHCONTENT_SOURCE_DIR_EIGEN ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
set( EIGEN_INCLUDE_INSTALL_DIR "include/feelpp/eigen" )

FetchContent_Declare(eigen GIT_REPOSITORY git@gitlab.com:feelpp/eigen.git GIT_TAG 3147391d946bb4b6c68edd901f2add6ac1f31f8c) # the commit hash for 1.5.0
FetchContent_MakeAvailable(eigen)
target_include_directories( feelpp_contrib BEFORE
      INTERFACE
      $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>
      $<BUILD_INTERFACE:${eigen_SOURCE_DIR}/unsupported>
      $<INSTALL_INTERFACE:include/feelpp/eigen>
      $<INSTALL_INTERFACE:include/feelpp/eigen/unsupported>
      )

target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_EIGEN3 )
target_link_libraries(feelpp_contrib INTERFACE Eigen3::Eigen )
# there is an issue with eigen+mkl
if ( MKL_FOUND )
  target_include_directories( eigen  INTERFACE ${MKL_INCLUDE_DIRS})
  target_compile_definitions(eigen INTERFACE EIGEN_USE_MKL_ALL)
  target_link_libraries(eigen INTERFACE ${MKL_LIBRARIES})  # error with cmake < 3.13, see CMP0079, use instead feelpp_contrib
  #target_link_libraries(feelpp_contrib INTERFACE ${MKL_LIBRARIES})
  #target_include_directories(eigen INTERFACE ${MKL_INCLUDE_DIRS})
  
endif()

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
message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring Feelpp Metis")
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



message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring NLopt")
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

message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring IPopt")
if( FEELPP_HAS_IPOPT )
  add_subdirectory(ipopt)
endif()


message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring PYBind11")
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
message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring Mongo")
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

message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring KWsys")
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

option( FEELPP_ENABLE_FMI4CPP "Enable FMI4CPP" ON )
if ( FEELPP_ENABLE_FMI4CPP )
  message(STATUS "------------------------------------------------------------")
  message(STATUS "Configuring FMI4CPP")
  set( USE_SYSTEM_CURL ON )
  set( BUILD_FMI4CPP_TESTS ON )
  set( FMI4CPP_HEADER_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/")
  set( FMI4CPP_CMAKE_INSTALL_DIR "share/feelpp/fmi4cpp")
  set( FMI4CPP_DOC_INSTALL_DIR "share/doc/feelpp/fmi4cpp")
  set(FETCHCONTENT_SOURCE_DIR_FMI4CPP ${CMAKE_CURRENT_SOURCE_DIR}/fmi4cpp)
  # patch includes fix for 18.04 compilation
  FetchContent_Declare(fmi4cpp GIT_REPOSITORY https://github.com/feelpp/fmi4cpp.git GIT_TAG 645b6a35b1c166ace666c5d387d5a88b074ace18) 

  FetchContent_MakeAvailable(fmi4cpp)
  add_library(feelpp_fmi4cpp INTERFACE)  
  set_target_properties(fmi4cpp PROPERTIES OUTPUT_NAME "feelpp_fmi4cpp")
  target_link_libraries(feelpp_fmi4cpp INTERFACE fmi4cpp)
  target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_FMI4CPP )
  add_library(Feelpp::feelpp_fmi4cpp ALIAS feelpp_fmi4cpp)
  target_include_directories(feelpp_fmi4cpp INTERFACE
      $<BUILD_INTERFACE:${FMI4CPP_SOURCE_DIR}/include/>
      $<INSTALL_INTERFACE:include/feelpp> )
  install(TARGETS feelpp_fmi4cpp fmi4cpp EXPORT feelpp-contrib-export-targets
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      )
  target_link_libraries( feelpp_contrib INTERFACE feelpp_fmi4cpp )
  set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} fmi4cpp/Contrib")
endif()

message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring FMILib")
option( FEELPP_ENABLE_FMILIB "Enable FMLIB C API" OFF )
# FMILib configure
# install rules.
set(FMILIB_BUILD_STATIC_LIB OFF)
set(FMILIB_BUILD_SHARED_LIB ON)
set(FMILIB_BUILD_TESTS ON)
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
  message(STATUS "------------------------------------------------------------")
  message(STATUS "Configuring GMSH")
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

message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring Tabulate")
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
message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring Indicators")
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
  
message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring Boost")
add_subdirectory( boost )
#include(feelpp.module.mongocxx)
#include(feelpp.module.hpddm)
#include(feelpp.module.nlopt)
#include(feelpp.module.ipopt)
#include(feelpp.module.cereal)
#include(feelpp.module.paralution)
#include(feelpp.module.jsonlab)

option( FEELPP_ENABLE_SIMPLE_WEB_SERVER "Enable Simple Web Server" OFF )
if ( FEELPP_ENABLE_SIMPLE_WEB_SERVER )
  message(STATUS "------------------------------------------------------------")
  message(STATUS "Configuring Simplee Web Server")
  feelppContribPrepare( simple-web-server )
  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    set(BUILD_TESTING OFF)
    add_subdirectory( simple-web-server )
    target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_SIMPLE_WEB_SERVER )
    target_link_libraries(feelpp_contrib INTERFACE simple-web-server )
    target_include_directories(feelpp_contrib INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/simple-web-server/>
      $<INSTALL_INTERFACE:include/feelpp> )
    set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} SimpleWebServer/Contrib")
  endif()
endif( FEELPP_ENABLE_SIMPLE_WEB_SERVER )

option( FEELPP_ENABLE_RANGE_V3 "Enable Range V3" ON )
if ( FEELPP_ENABLE_RANGE_V3 )
  message(STATUS "------------------------------------------------------------")
  message(STATUS "Configuring Range V3")
  feelppContribPrepare( range-v3 )
  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    set( RANGES_CXX_STD "17" )
    if ( NOT ${CMAKE_INSTALL_INCLUDEDIR} MATCHES "feelpp$")
      set( CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/feelpp")
    endif()
    message( STATUS "range-v3 CMAKE_INSTALL°INCLUDE_DIR: ${CMAKE_INSTALL_INCLUDEDIR}")
    add_subdirectory( range-v3 )
    target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_RANGE_V3 )
    target_link_libraries(feelpp_contrib INTERFACE range-v3 )
    target_include_directories(feelpp_contrib INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/range-v3/>
      $<INSTALL_INTERFACE:include/feelpp> )
    set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} range-v3/Contrib")
  endif()
endif( FEELPP_ENABLE_RANGE_V3 )

option( FEELPP_ENABLE_MATPLOT "Enable MatPlot++" OFF )
if ( FEELPP_ENABLE_MATPLOT )
  message(STATUS "------------------------------------------------------------")
  message(STATUS "Configuring MatPlot++")
  feelppContribPrepare( matplot )
  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    if ( NOT ${CMAKE_INSTALL_INCLUDEDIR} MATCHES "feelpp$")
      set( CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/feelpp")
    endif()
    message( STATUS "matplot++ CMAKE_INSTALL_INCLUDE_DIR: ${CMAKE_INSTALL_INCLUDEDIR}")
    set( BUILD_SHARED_LIBS ON )
    set( BUILD_INSTALLER ON )
    set(OPTIMIZATION_FLAG_IS_SET ON)
    add_subdirectory( matplot )
    target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_MATPLOT )
    set_target_properties(matplot PROPERTIES OUTPUT_NAME "feelpp_matplot")
    target_link_libraries( feelpp_contrib INTERFACE matplot )
    install(TARGETS matplot EXPORT feelpp-contrib-export-targets
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      )
#    target_link_libraries(feelpp_contrib INTERFACE matplot )
    target_include_directories(feelpp_contrib INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/contrib/matplot/source/>
      $<INSTALL_INTERFACE:include/feelpp> )
    set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} MatPPlot++/Contrib")
  endif()
endif( FEELPP_ENABLE_MATPLOT )

option( FEELPP_PARMMG_DOWNLOAD_METIS "Download Metis for PARMMG" OFF )
option( FEELPP_ENABLE_MMG "Enable MMG/PARMMG" ON )
option( FEELPP_PARMMG_DOWNLOAD_METIS "Disable metis dowload" OFF )
if ( FEELPP_ENABLE_MMG )
message(STATUS "------------------------------------------------------------")
message(STATUS "Configuring Mmg")
feelppContribPrepare( mmg )
if( FEELPP_CONTRIB_PREPARE_SUCCEED )
  set(MMG_CI OFF)
  set(MMG2D_CI OFF)
  set(MMG3D_CI OFF)
  set(MMGS_CI OFF)
  set(BUILD_SHARED_LIBS ON)
  set(ONLY_VERY_SHORT_TESTS ON)
  set(MMG_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/feel/cmake/modules")
  if ( NOT ${CMAKE_INSTALL_INCLUDEDIR} MATCHES "feelpp$")
    set(CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/feelpp" )
  endif()
  message( STATUS "mmg CMAKE_INSTALL°INCLUDE_DIR: ${CMAKE_INSTALL_INCLUDEDIR}")
  add_subdirectory( mmg )
  target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_MMG )
  target_link_libraries(feelpp_contrib INTERFACE libfeelpp_mmg_so )
  set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Mmg/Contrib")
  if ( SCOTCH_FOUND )
    target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_SCOTCH )
  endif() 
  

  # setup parmmg only if mmg is ok
  feelppContribPrepare( parmmg )
  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    SET ( LIBPARMMG_STATIC OFF )
    SET ( LIBPARMMG_SHARED ON )
    SET( DOWNLOAD_MMG OFF )
    SET( DOWNLOAD_METIS ${FEELPP_PARMMG_DOWNLOAD_METIS} )
    
#    SET(MMG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mmg)
#    SET(MMG_BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/mmg)
#    SET(MMG_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/mmg)
    add_subdirectory( parmmg )
    target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_PARMMG )
    target_link_libraries(feelpp_contrib INTERFACE libfeelpp_parmmg_so )
    set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} ParMmg/Contrib")

  endif()
endif(FEELPP_CONTRIB_PREPARE_SUCCEED)
endif(FEELPP_ENABLE_MMG)

option( FEELPP_ENABLE_CPR "Enable CPR" ON )
if ( FEELPP_ENABLE_CPR )
  message(STATUS "------------------------------------------------------------")
  message(STATUS "Configuring CPR")
  set(FETCHCONTENT_SOURCE_DIR_CPR ${CMAKE_CURRENT_SOURCE_DIR}/cpr)
  # patch includes fix for 18.04 compilation
  FetchContent_Declare(cpr GIT_REPOSITORY https://github.com/feelpp/cpr.git GIT_TAG e554da6bde33edc35f1c28bd8197ae50bdf4a4d4) 
  set( USE_SYSTEM_CURL ON )
  set( BUILD_CPR_TESTS OFF )
  if ( NOT ${CMAKE_INSTALL_INCLUDEDIR} MATCHES "feelpp$")
      set( CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/feelpp")
  endif()
  FetchContent_MakeAvailable(cpr)
  add_library(feelpp_cpr INTERFACE)  
  set_target_properties(cpr PROPERTIES OUTPUT_NAME "feelpp_cpr")
  target_link_libraries(feelpp_cpr INTERFACE cpr)
  target_compile_definitions( feelpp_contrib INTERFACE FEELPP_HAS_CPR )
  add_library(Feelpp::feelpp_cpr ALIAS feelpp_cpr)
  target_include_directories(feelpp_cpr INTERFACE
      $<BUILD_INTERFACE:${cpr_SOURCE_DIR}/include/>
      $<INSTALL_INTERFACE:include/feelpp> )
  install(TARGETS feelpp_cpr cpr EXPORT feelpp-contrib-export-targets
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      )
  target_link_libraries( feelpp_contrib INTERFACE feelpp_cpr )
  set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} CPR/Contrib")
endif()

# ----------------------------------------------------------
# A fast and simple to use logging library.
set(FETCHCONTENT_SOURCE_DIR_SPDLOG ${CMAKE_CURRENT_SOURCE_DIR}/spdlog)
FetchContent_Declare(   
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.8.2)

FetchContent_GetProperties(spdlog)
if(NOT spdlog_POPULATED)
    message(STATUS "------------------------------------------------------------")
    message(STATUS "Configuring SPDLOG")
    set(SPDLOG_FMT_EXTERNAL ON)
    #FetchContent_Populate(spdlog)
    FetchContent_MakeAvailable(spdlog)
#    add_subdirectory(
#    ${spdlog_SOURCE_DIR}
#    ${spdlog_BINARY_DIR})
    target_link_libraries( feelpp_contrib INTERFACE spdlog::spdlog_header_only )
    install(TARGETS spdlog_header_only EXPORT feelpp-contrib-export-targets
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      )
endif() 

option( FEELPP_ENABLE_TASKFLOW "Enable TASKFLOW" ON )
if ( FEELPP_ENABLE_TASKFLOW )
  message(STATUS "------------------------------------------------------------")
  message(STATUS "Configuring TASKFLOW")
#  set(FETCHCONTENT_SOURCE_DIR_TASKFLOW ${CMAKE_CURRENT_SOURCE_DIR}/taskflow)
  # patch includes fix for 18.04 compilation
  set(FETCHCONTENT_SOURCE_DIR_TASKFLOW ${CMAKE_CURRENT_SOURCE_DIR}/taskflow)
  FetchContent_Declare(taskflow GIT_REPOSITORY https://github.com/feelpp/taskflow.git GIT_TAG 31699d74dfb92d4d55fc0f3c87899c6ad8b51d80)
  set( BUILD_TASKFLOW_TESTS OFF )
  set(TF_CMAKE_DIR  "${CMAKE_INSTALL_DATADIR}/feel/cmake/modules")
  FetchContent_MakeAvailable(taskflow)
  add_library(feelpp_taskflow INTERFACE)
  set_target_properties(Taskflow PROPERTIES OUTPUT_NAME "feelpp_taskflow")
  target_link_libraries(feelpp_taskflow INTERFACE Taskflow spdlog::spdlog_header_only)
  add_library(Feelpp::feelpp_taskflow ALIAS feelpp_taskflow)
  target_include_directories(feelpp_taskflow INTERFACE
      $<BUILD_INTERFACE:${taskflow_SOURCE_DIR}/>
      $<INSTALL_INTERFACE:include/feelpp> )
  install(TARGETS feelpp_taskflow Taskflow EXPORT feelpp-contrib-export-targets
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      )
  target_link_libraries( feelpp_contrib INTERFACE feelpp_taskflow )
  set(FEELPP_ENABLED_OPTIONS_CONTRIB "${FEELPP_ENABLED_OPTIONS_CONTRIB} Taskflow")
endif()




# Add an info message to be displayed at the end of the cmake process.
if( FEELPP_CONTRIB_SUBMODULE_UPDATED )
  list( APPEND FEELPP_MESSAGE_INFO_END "Feel++ submodules already initialized!\nPlease make sure submodules are up to date (run `git submodule update --init --recursive` in source directory)" )
  set( FEELPP_MESSAGE_INFO_END ${FEELPP_MESSAGE_INFO_END} )
endif()


INSTALL(TARGETS feelpp_contrib DESTINATION lib/ COMPONENT Libs EXPORT feelpp-contrib-export-targets)
install(EXPORT feelpp-contrib-export-targets NAMESPACE Feelpp:: DESTINATION share/feelpp/feel/cmake/modules)

message(STATUS "Contrib libraries configured")
message(STATUS "============================================================")
