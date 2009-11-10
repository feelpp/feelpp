cmake_minimum_required (VERSION 2.6)

SET( CMAKE_CXX_FLAGS "-pipe -Wall -O2 ")
SET( CMAKE_C_FLAGS "-pipe -Wall -O2")

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
set(LIFE_VERSION_MICRO "11")
set(LIFE_REVISION "")
set(LIFE_BUILDID "")
set(LIFE_VERSION "(((${LIFE_VERSION_MAJOR}) << 16) | ((${LIFE_VERSION_MINOR}) << 9) | (${LIFE_VERSION_MICRO}))" )


SET(LIFE_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(LIFE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")

SET( CMAKE_MODULE_PATH
  ${LIFE_HOME_DIR}/cmake
  #    ${LIFE_HOME_DIR}/cmake/utils
  #    ${LIFE_HOME_DIR}/cmake/package_arch
  #    ${LIFE_HOME_DIR}/cmake/config_tests
  )

OPTION(ENABLE_INSTANTIATION_MODE "Instantiation mode" ON )
MARK_AS_ADVANCED(ENABLE_INSTANTIATION_MODE)
IF ( ENABLE_INSTANTIATION_MODE )
ADD_DEFINITIONS( -DLIFE_INSTANTIATION_MODE=1 )
ENDIF()

set(LIFE_MESH_MAX_ORDER "5")

#INCLUDE(PackageArchGlobalMacros)
# INCLUDE(LifeGlobalMacros)
#INCLUDE(AdvancedSet)
#INCLUDE(AdvancedOption)


INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckLibraryExists)


Include(FindMPI)
FIND_PACKAGE(MPI)
SET(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")
IF ( MPI_FOUND )
  ADD_DEFINITIONS( -DHAVE_MPI -DHAVE_MPI_H )
ENDIF()

#
# Xml2
#
Include(FindLibXml2)
FIND_PACKAGE(LibXml2)
#
# Blas and Lapack
#
INCLUDE(FindLAPACK)
FIND_PACKAGE(LAPACK)


Include(FindBoost)
FIND_PACKAGE(Boost COMPONENTS date_time filesystem system program_options unit_test_framework signals  mpi regex  serialization)
SET(Boost_ADDITIONAL_VERSIONS "1.39" )
ADD_DEFINITIONS( -DBOOST_PARAMETER_MAX_ARITY=15 -DBOOST_TEST_DYN_LINK )

#
# Scotch
#
#CheckIncludeFileCXX( ptscotch.h HAVE_PTSCOTCH_H )
#CheckIncludeFileCXX( scotch.h HAVE_SCOTCH_H )

#
# Petsc
#
INCLUDE(FindPETSc)
find_package( PETSc )
if ( PETSC_FOUND )
  add_definitions( -DHAVE_PETSC -DHAVE_PETSC_H )
endif()
FIND_PATH(
  SLEPC_INCLUDE_DIR
  NAMES slepc.h
  PATHS /usr/include/slepc
  )
SET(CMAKE_REQUIRED_INCLUDES "${SLEPC_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")
CHECK_INCLUDE_FILE( slepc.h HAVE_SLEPC_H )
FIND_LIBRARY(SLEPC_LIB_SLEPC     slepc )
SET(SLEPC_LIBRARIES ${SLEPC_LIB_SLEPC})



#
# Trilinos
#
OPTION(ENABLE_TRILINOS "enable trilinos" ON)
#IF ( ENABLE_TRILINOS )
INCLUDE(FindTrilinos)
INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIR})
#ENDIF( ENABLE_TRILINOS )

INCLUDE(FindVTK)
FIND_PACKAGE(VTK)

# generate configuration header
CONFIGURE_FILE(lifeconfig.h.in lifeconfig.h)

#
#
#
INCLUDE_DIRECTORIES (
  ${LIFE_SOURCE_DIR}/
  ${LIFE_SOURCE_DIR}/contrib/gmm/include

  ${SLEPC_INCLUDE_DIR}
  ${PETSC_INCLUDE_DIR}
  ${PETSCCONF_INCLUDE_DIR}

  ${MPI_INCLUDE_PATH}
  ${BOOST_INCLUDE_PATH}
  ${LIBXML2_INCLUDE_DIR}
  )

LINK_DIRECTORIES(
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


ADD_SUBDIRECTORY ( life )

SET(LIFE_LIBRARIES
  life
  lifematerial
  lifefilters
  lifediscr
  lifealg
  lifemesh
  lifecore
  ${TRILINOS_LIBRARIES}
  ${PETSC_LIBRARIES}
  ${Boost_LIBRARIES}
  ${MPI_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LIBXML2_LIBRARIES}
  )

ADD_SUBDIRECTORY ( doc )
ADD_SUBDIRECTORY ( examples )

#
# Enable testing
#
macro(life_add_test testname)

  set(targetname test_${testname})
  set(filename test_${testname}.cpp)
  add_executable(${targetname} ${filename})
  target_link_libraries(${targetname} ${LIFE_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
  add_test( ${testname} "${targetname}" )
endmacro(life_add_test)



ENABLE_TESTING()
foreach( testdir core alg mesh poly filters material discr )
  add_subdirectory( testsuite/life${testdir} )
endforeach()

#
# Packaging
#
include( CPack )
