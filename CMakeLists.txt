cmake_minimum_required (VERSION 2.6)

SET(CMAKE_CXX_FLAGS "-std=c++0x" CACHE STRING "Default flags" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "-std=c++0x ${CMAKE_CXX_FLAGS_DEBUG} " CACHE STRING "Debug flags" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x ${CMAKE_CXX_FLAGS_RELEASE} " CACHE STRING "Release flags" FORCE)

#SET( CMAKE_CXX_FLAGS "-pipe -Wall -O2 ")
#SET( CMAKE_C_FLAGS "-pipe -Wall -O2")

project (Life C CXX Fortran)

# FIND_PACKAGE(Subversion)
# IF(Subversion_FOUND)
#   Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
#   MESSAGE("Current revision is ${Project_WC_REVISION}")
#   Subversion_WC_LOG(${PROJECT_SOURCE_DIR} Project)
#   MESSAGE("Last changed log is ${Project_LAST_CHANGED_LOG}")
# ENDIF(Subversion_FOUND)

set(LIFE_VERSION_MAJOR "0")
set(LIFE_VERSION_MINOR "9")
set(LIFE_VERSION_MICRO "24")
set(LIFE_VERSION_EXTRA "0")
set(LIFE_REVISION "1")
set(LIFE_BUILDID "1")
set(LIFE_VERSION "(((${LIFE_VERSION_MAJOR}) << 16) | ((${LIFE_VERSION_MINOR}) << 9) | (${LIFE_VERSION_MICRO}))" )
set(LIFE_VERSION_STRING "${LIFE_VERSION_MAJOR}.${LIFE_VERSION_MINOR}.${LIFE_VERSION_MICRO}")
set(LIFE_SHARED_VERSION "1.0.0" )
set(LIFE_SHARED_SOVERSION "1" )

SET(LIFE_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(LIFE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(LIFE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")

SET( CMAKE_MODULE_PATH
  ${LIFE_HOME_DIR}/cmake
  #    ${LIFE_HOME_DIR}/cmake/utils
  #    ${LIFE_HOME_DIR}/cmake/package_arch
  #    ${LIFE_HOME_DIR}/cmake/config_tests
  )

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

OPTION(ENABLE_MOVE_SEMANTICS "enable move semantics(elision)" ON )
MARK_AS_ADVANCED(ENABLE_MOVE_SEMANTICS)
IF ( ENABLE_MOVE_SEMANTICS )
  SET( BOOST_UBLAS_MOVE_SEMANTICS 1 CACHE STRING "Enable Boost Ublas move semantics" FORCE )
  ADD_DEFINITIONS( -DBOOST_UBLAS_MOVE_SEMANTICS )
ENDIF( ENABLE_MOVE_SEMANTICS)

OPTION(ENABLE_INSTANTIATION_MODE "Instantiation mode" ON )
MARK_AS_ADVANCED(ENABLE_INSTANTIATION_MODE)
IF ( ENABLE_INSTANTIATION_MODE )
  SET( LIFE_INSTANTIATION_MODE 1 )
ENDIF()

SET(LIFE_MESH_MAX_ORDER "5")

#INCLUDE(PackageArchGlobalMacros)
# INCLUDE(LifeGlobalMacros)
#INCLUDE(AdvancedSet)
#INCLUDE(AdvancedOption)


INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckLibraryExists)
# INCLUDE(parse_arguments)

# MACRO(CAR var)
#   SET(${var} ${ARGV1})
# ENDMACRO(CAR)

# MACRO(CDR var junk)
#   SET(${var} ${ARGN})
# ENDMACRO(CDR)
INCLUDE(CheckTypeSize)
CHECK_TYPE_SIZE(int SIZE_INT )
CHECK_TYPE_SIZE(uint SIZE_UINT )
CHECK_TYPE_SIZE(size_t SIZE_SIZE_T )
CHECK_TYPE_SIZE(long SIZE_LONG )
CHECK_TYPE_SIZE(float SIZE_FLOAT )
CHECK_TYPE_SIZE(double SIZE_DOUBLE )
CHECK_TYPE_SIZE("long double" SIZE_LONG_DOUBLE)
MESSAGE(STATUS "SIZE_INT=${SIZE_INT}")
MESSAGE(STATUS "SIZE_UINT=${SIZE_UINT}")
MESSAGE(STATUS "SIZE_SIZE_T=${SIZE_SIZE_T}")
MESSAGE(STATUS "SIZE_LONG=${SIZE_LONG}")
MESSAGE(STATUS "SIZE_FLOAT=${SIZE_FLOAT}")
MESSAGE(STATUS "SIZE_DOUBLE=${SIZE_DOUBLE}")
MESSAGE(STATUS "SIZE_LONG_DOUBLE=${SIZE_LONG_DOUBLE}")

FIND_PACKAGE(OpenMP)

FIND_PACKAGE(MPI)
IF ( MPI_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")
  SET( HAVE_MPI 1 )
  SET( HAVE_MPI_H 1 )
  ADD_DEFINITIONS( -DHAVE_MPI=1 -DHAVE_MPI_H=1 )
  SET(LIFE_BOOST_MPI mpi)
ENDIF()



Check_Include_File_CXX(dlfcn.h HAVE_DLFCN_H)
if ( HAVE_DLFCN_H )
  add_definitions(-DHAVE_DLFCN_H)
endif()
CHECK_LIBRARY_EXISTS (dl dlopen "" HAVE_LIBDL)
IF (HAVE_LIBDL)
  set(DL_LIBRARY dl)
  SET (HAVE_DLOPEN 1)
  add_definitions(-DHAVE_DLOPEN)
ELSE ()
  CHECK_FUNCTION_EXISTS (dlopen HAVE_DLOPEN)
ENDIF (HAVE_LIBDL)

#
# Xml2
#
FIND_PACKAGE(LibXml2 REQUIRED)

#
# Blas and Lapack
#
FIND_PACKAGE(LAPACK)

FIND_PACKAGE(Boost COMPONENTS date_time filesystem system program_options unit_test_framework signals  ${LIFE_BOOST_MPI} regex  serialization)
set(Boost_ADDITIONAL_VERSIONS "1.39" "1.40" "1.41" "1.42" )
set( BOOST_PARAMETER_MAX_ARITY 15 )
add_definitions( -DBOOST_PARAMETER_MAX_ARITY=${BOOST_PARAMETER_MAX_ARITY} -DBOOST_TEST_DYN_LINK )

#
# Scotch
#
#CheckIncludeFileCXX( ptscotch.h HAVE_PTSCOTCH_H )
#CheckIncludeFileCXX( scotch.h HAVE_SCOTCH_H )
#
# Metis
#
FIND_PACKAGE(Metis)
if ( METIS_FOUND )
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
  LINK_DIRECTORIES(${METIS_LIBRARIES})
endif( METIS_FOUND )

#
# Petsc
#
FIND_PACKAGE( PETSc )
if ( PETSC_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDES};${CMAKE_REQUIRED_INCLUDES}")
endif( PETSC_FOUND )

#
# SLEPc
#
FIND_PACKAGE( SLEPc )
if ( SLEPC_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${SLEPC_INCLUDES};${CMAKE_REQUIRED_INCLUDES}")
  INCLUDE_DIRECTORIES( ${SLEPC_INCLUDE_DIR} )
endif(SLEPC_FOUND)

#
# parpack
#
FIND_LIBRARY(PARPACK_LIBRARY NAMES parpack)
if (PARPACK_LIBRARY)
    SET(PARPACK_LIBRARIES ${PARPACK_LIBRARY})
endif()
MARK_AS_ADVANCED( PARPACK_LIBRARY )


#
# Trilinos
#
FIND_PACKAGE(Trilinos)
if ( TRILINOS_FOUND )
  INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIR})
endif( TRILINOS_FOUND )

#
# OpenTURNS
#
FIND_PACKAGE( OpenTURNS )
if ( OpenTURNS_FOUND )
  INCLUDE_DIRECTORIES(${OpenTURNS_INCLUDE_DIRS})
endif( OpenTURNS_FOUND )

#
# VTK
#
FIND_PACKAGE(VTK)
if ( VTK_FOUND )
  set(HAVE_VTK 1)
  SET(VTK_LIBRARIES "-lvtkRendering -lvtkGraphics -lvtkImaging -lvtkCommon" )
  INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
  MARK_AS_ADVANCED( VTK_DIR )
endif()

# generate configuration header
CONFIGURE_FILE(lifeconfig.h.in lifeconfig.h  )

#
#
#
INCLUDE_DIRECTORIES (
  ${LIFE_BUILD_DIR}/
  ${LIFE_SOURCE_DIR}/
  ${LIFE_SOURCE_DIR}/contrib/gmm/include

  ${PETSC_INCLUDE_DIR}
  ${PETSCCONF_INCLUDE_DIR}

  ${MPI_INCLUDE_PATH}
  ${BOOST_INCLUDE_PATH}
  ${LIBXML2_INCLUDE_DIR}
  )

LINK_DIRECTORIES(
  ${VTK_LIBRARY_DIRS}
  ${BOOST_LIBRARY_PATH}
  ${MPI_LIBRARY_PATH}
  )

#
# programs
#
find_program( GMSH gmsh DOC "Gmsh mesh generator" )
if ( GMSH )
  ADD_DEFINITIONS( -DHAVE_GMSH=1 )
endif()
mark_as_advanced( GMSH )

# ${LIFE_SOURCE_DIR}/life/lifecore
# ${LIFE_SOURCE_DIR}/life/lifealg
# ${LIFE_SOURCE_DIR}/life/lifemesh
# ${LIFE_SOURCE_DIR}/life/lifepoly
# ${LIFE_SOURCE_DIR}/life/lifediscr
# ${LIFE_SOURCE_DIR}/life/lifevf
# ${LIFE_SOURCE_DIR}/life/lifefilters
# ${LIFE_SOURCE_DIR}/life/lifesystem

#
# Python
#
FIND_PACKAGE(PythonInterp 2.2 REQUIRED)

#
# Life
#
ADD_SUBDIRECTORY ( life )

SET(LIFE_LIBRARIES
  life
  lifematerial
  lifefilters
  lifediscr
  lifealg
  lifemesh
  lifecore
  ${OT_LIBRARIES}
  ${TRILINOS_LIBRARIES}
  ${VTK_LIBRARIES}
  ${SLEPC_LIBRARIES}
  ${PETSC_LIBRARIES}
  ${Boost_LIBRARIES}
  ${METIS_LIBRARY}
  ${MPI_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LIBXML2_LIBRARIES}
  ${PARPACK_LIBRARIES}
  ${DL_LIBRARY}
  )



OPTION(LIFE_ENABLE_DOCUMENTATION "enable life documentation" OFF)
OPTION(LIFE_ENABLE_BENCHMARKS "enable life benchmarks" OFF)
OPTION(LIFE_ENABLE_EXAMPLES "enable life examples" OFF)
OPTION(LIFE_ENABLE_APPLICATIONS "enable life applications" OFF)
OPTION(LIFE_ENABLE_TESTS "enable life tests" OFF)
OPTION(LIFE_ENABLE_ALL "enable all life modules" OFF)

if ( LIFE_ENABLE_ALL )
  set( LIFE_ENABLE_DOCUMENTATION ON )
  set( LIFE_ENABLE_BENCHMARKS ON )
  set( LIFE_ENABLE_EXAMPLES ON )
  set( LIFE_ENABLE_APPLICATIONS ON )
  set( LIFE_ENABLE_TESTS ON )
endif( LIFE_ENABLE_ALL )

#
# Enable testing
#
INCLUDE(CTest)
ENABLE_TESTING()

if ( LIFE_ENABLE_DOCUMENTATION )
  ADD_SUBDIRECTORY ( doc )
endif ( LIFE_ENABLE_DOCUMENTATION )

if ( BUILD_TESTING OR LIFE_ENABLE_BENCHMARKS )
  ADD_SUBDIRECTORY ( benchmarks )
endif()


if ( LIFE_ENABLE_EXAMPLES )
 ADD_SUBDIRECTORY ( examples )
endif()

if ( LIFE_ENABLE_APPLICATIONS )
 ADD_SUBDIRECTORY ( applications )
 if ( EXISTS ${LIFE_SOURCE_DIR}/research )
  ADD_SUBDIRECTORY ( research )
 endif()
endif( LIFE_ENABLE_APPLICATIONS )




IF(BUILD_TESTING OR LIFE_ENABLE_TESTS )

  macro(life_add_test testname)
    set(targetname test_${testname})
    set(filename test_${testname}.cpp)
    add_executable(${targetname} ${filename})
    target_link_libraries(${targetname} ${LIFE_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
    add_test( test_${testname} "${targetname}" )
  endmacro(life_add_test)

  add_subdirectory( testsuite )

endif()
#
# Packaging
#
INCLUDE(InstallRequiredSystemLibraries)


foreach(includedir lifecore lifealg lifemesh lifepoly lifefilters lifediscr lifevf lifematerial lifesystem )
  FILE(GLOB files "life/${includedir}/*.hpp")
  if ( NOT ENABLE_INSTANTIATION_MODE )
    FILE(GLOB cppfiles "life/${includedir}/*.cpp")
  endif()
  INSTALL(FILES ${files} ${cppfiles} DESTINATION include/life/${includedir})
endforeach()
FILE(GLOB files "life/*.hpp")
INSTALL(FILES ${files} DESTINATION include/life)
FILE(GLOB files "contrib/gmm/include/*.h")
INSTALL(FILES ${files} DESTINATION include/life)
FILE(GLOB files ${CMAKE_CURRENT_BINARY_DIR}/lifeconfig.h)
INSTALL(FILES ${files} DESTINATION include/)
foreach(libdir lifecore lifealg lifemesh lifefilters lifevf
lifematerial lifediscr)
  FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/life/${libdir}/*.so*")
  INSTALL(FILES ${files} DESTINATION lib/)
endforeach()
FILE(GLOB libso "${CMAKE_CURRENT_BINARY_DIR}/life/*.so*")
INSTALL(FILES ${libso} DESTINATION lib/)
foreach( programs_dir examples/generic examples/heat examples/solid doc/tutorial applications/polyvis benchmarks/convergence )
  FILE(GLOB programs "${CMAKE_CURRENT_BINARY_DIR}/${programs_dir}/life_*")
  INSTALL(FILES ${programs} DESTINATION bin/)
endforeach()
# documentation and examples
FILE(GLOB pdfs "${CMAKE_CURRENT_BINARY_DIR}/doc/tutorial/life-manual.pdf")
INSTALL(FILES ${pdfs} DESTINATION share/doc/life/)
FILE(GLOB examples "${CMAKE_CURRENT_SOURCE_DIR}/doc/tutorial/*.*pp")
INSTALL(FILES ${examples} DESTINATION share/doc/life/examples/)


##
## Archive generation using cpack
##
SET(CPACK_GENERATOR "TGZ")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "life")
SET(CPACK_PACKAGE_VENDOR "Christophe Prud'homme")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "${LIFE_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${LIFE_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${LIFE_VERSION_MICRO}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "life")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "CPackSourceConfig.cmake")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "life-${LIFE_VERSION_MAJOR}.${LIFE_VERSION_MINOR}.${LIFE_VERSION_MICRO}")
SET(CPACK_SOURCE_STRIP_FILES "")
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
set(CPACK_SOURCE_IGNORE_FILES
  "/\\\\.git/;\\\\.gitignore;/\\\\.svn;"
  "/admin/;/tools/;/Templates/;"
  "/auto/;/auto.*/;"
  "/TAGS;/#.*;/.*~;/.cvsignore;/.bzrignore"
  "${PROJECT_SOURCE_DIR}/benchmarks/turek/"
  "${PROJECT_SOURCE_DIR}/benchmarks/stokes/"
  "${PROJECT_SOURCE_DIR}/benchmarks/ethiersteinman/"
  "${PROJECT_SOURCE_DIR}/benchmarks/kovasznay/"
  "${PROJECT_SOURCE_DIR}/doc/poster/"
  "${PROJECT_SOURCE_DIR}/doc/tutorial/pdfs/"
  "${PROJECT_SOURCE_DIR}/doc/figures/backgrounds/"
  "${PROJECT_SOURCE_DIR}/doc/figures/logos/"
  "${PROJECT_SOURCE_DIR}/examples/fluid/"
  "${PROJECT_SOURCE_DIR}/examples/levelset/"
  "${PROJECT_SOURCE_DIR}/examples/pbeq/"
  "${PROJECT_SOURCE_DIR}/applications/opus/"
  "${PROJECT_SOURCE_DIR}/research/"  )

include( CPack )


