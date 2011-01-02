# - Find Feel
# This module looks for Feel (Library for the Finite Element Method) support
# it will define the following values
#  FEEL_INCLUDE_DIR = where feel/feelcore/feel.hpp can be found
#  FEEL_LIBRARY    = the library to link in

INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckLibraryExists)


FIND_PACKAGE(MPI)
SET(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")
IF ( MPI_FOUND )
  ADD_DEFINITIONS( -DHAVE_MPI -DHAVE_MPI_H )
ENDIF()

#
# Xml2
#
FIND_PACKAGE(LibXml2)
#
# Blas and Lapack
#
FIND_PACKAGE(LAPACK)


find_package(Boost COMPONENTS date_time filesystem system program_options unit_test_framework signals  mpi regex  serialization)
set(Boost_ADDITIONAL_VERSIONS "1.39" "1.40" "1.41")
set( BOOST_PARAMETER_MAX_ARITY 15 )
add_definitions( -DBOOST_PARAMETER_MAX_ARITY=${BOOST_PARAMETER_MAX_ARITY} -DBOOST_TEST_DYN_LINK )

#
# Scotch
#
#CheckIncludeFileCXX( ptscotch.h HAVE_PTSCOTCH_H )
#CheckIncludeFileCXX( scotch.h HAVE_SCOTCH_H )
# Metis
#
FIND_PACKAGE(Metis)
if ( METIS_FOUND )
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
  LINK_DIRECTORIES(${METIS_LIBRARIES})
endif()

#
# Petsc
#
FIND_PACKAGE( PETSc REQUIRED )
if ( PETSC_FOUND )
  add_definitions( -DHAVE_PETSC -DHAVE_PETSC_H )
  MARK_AS_ADVANCED( PETSC_CURRENT PETSC_DIR PETSC_ARCH )
endif()

find_path (SLEPC_DIR include/slepc.h
  HINTS ENV SLEPC_DIR
  PATHS
  /usr/lib/slepcdir/3.0.0 # Debian
  $ENV{HOME}/slepc
  DOC "SLEPc Directory")
SET(CMAKE_REQUIRED_INCLUDES "${SLEPC_DIR}/include;${CMAKE_REQUIRED_INCLUDES}")
CHECK_INCLUDE_FILE( slepc.h HAVE_SLEPC_H )
FIND_LIBRARY(SLEPC_LIB_SLEPC     slepc )
SET(SLEPC_LIBRARIES ${SLEPC_LIB_SLEPC})
# message( "*** SLEPc directory : ${SLEPC_DIR}" )
if (HAVE_SLEPC_H AND SLEPC_DIR)
  set(HAVE_SLEPC 1)
  MARK_AS_ADVANCED( SLEPC_DIR SLEPC_LIB_SLEPC )
endif()


#
# parpack
#
FIND_LIBRARY(PARPACK_LIBRARY NAMES parpack)
if (PARPACK_LIBRARY)
  SET(PARPACK_LIBRARIES ${PARPACK_LIBRARY})
  MARK_AS_ADVANCED( PARPACK_LIBRARY )
endif()



#
# Trilinos
#
find_package(Trilinos)
if ( TRILINOS_FOUND )
  INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIR})
endif()

#
# OpenTURNS
#
FIND_PACKAGE( OpenTURNS )
#if ( OpenTURNS_FOUND )
INCLUDE_DIRECTORIES(${OpenTURNS_INCLUDE_DIRS})
ADD_DEFINITIONS( ${OpenTURNS_WRAPPER_DEFINITIONS})
#endif( OpenTURNS_FOUND )

#
# VTK
#
find_package(VTK)
if ( VTK_FOUND )
  set(HAVE_VTK 1)
  SET(VTK_LIBRARIES "-lvtkRendering -lvtkGraphics -lvtkImaging -lvtkCommon" )
  MARK_AS_ADVANCED( VTK_DIR )
endif()

#
# gmsh
#
find_program( GMSH gmsh DOC "Gmsh mesh generator" )
if ( GMSH )
  ADD_DEFINITIONS( -DHAVE_GMSH=1 )
endif()
mark_as_advanced( GMSH )

#
# Feel
#
SET(FEEL_FOUND 0)
FIND_PATH(FEEL_INCLUDE_DIR feelconfig.h  PATHS /usr/include/feel /usr/lib/feel/include /opt/feel/include /usr/ljk/include/feel /usr/local  )

ADD_DEFINITIONS( -DFEEL_INSTANTIATION_MODE=1 )

FIND_LIBRARY(FEEL_LIBRARY        feel        PATHS /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
FIND_LIBRARY(FEELCORE_LIBRARY    feelcore    PATHS /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
FIND_LIBRARY(FEELALG_LIBRARY     feelalg     PATHS /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
FIND_LIBRARY(FEELMESH_LIBRARY    feelmesh    PATHS /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
FIND_LIBRARY(FEELDISCR_LIBRARY   feeldiscr   PATHS /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
FIND_LIBRARY(FEELFILTERS_LIBRARY feelfilters PATHS /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
FIND_LIBRARY(FEELMATERIAL_LIBRARY feelmaterial PATHS /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )

SET(FEEL_LIBRARIES
#  ${FEEL_LIBRARY}
#  ${FEELMATERIAL_LIBRARY}
#  ${FEELFILTERS_LIBRARY}
#  ${FEELDISCR_LIBRARY}
#  ${FEELALG_LIBRARY}
#  ${FEELMESH_LIBRARY}
#  ${FEELCORE_LIBRARY}
feel
feelmaterial
feelfilters
feeldiscr
feelalg
feelmesh
feelcore
  ${OT_LIBRARIES}
  ${TRILINOS_LIBRARIES}
  ${VTK_LIBRARIES}
  ${SLEPC_LIBRARIES}
  ${PETSC_LIBRARIES}
  ${Boost_LIBRARIES}
  ${METIS_LIBRARIES}
  ${MPI_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LIBXML2_LIBRARIES}
  ${PARPACK_LIBRARIES}
  )

INCLUDE_DIRECTORIES (
  ${FEEL_INCLUDE_DIR}

  ${OT_INCLUDE_DIR}

  ${SLEPC_INCLUDE_DIR}
  ${PETSC_INCLUDE_DIR}
  ${PETSCCONF_INCLUDE_DIR}

  ${VTK_INCLUDE_DIRS}

  ${MPI_INCLUDE_PATH}
  ${BOOST_INCLUDE_PATH}
  ${LIBXML2_INCLUDE_DIR}
  )
IF(FEEL_INCLUDE_DIR AND FEEL_LIBRARIES)
   SET(FEEL_FOUND 1)
ENDIF(FEEL_INCLUDE_DIR AND FEEL_LIBRARIES)

MARK_AS_ADVANCED(
  FEEL_LIBRARY
  FEELCORE_LIBRARY
  FEELALG_LIBRARY
  FEELMESH_LIBRARY
  FEELDISCR_LIBRARY
  FEELFILTERS_LIBRARY

  BOOSTMPI_LIBRARY
  BOOSTSERIALIZATION_LIBRARY
  )
