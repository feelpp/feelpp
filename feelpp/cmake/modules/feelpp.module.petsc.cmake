if ( APPLE )
  message(STATUS "clang rt: ${PETSC_CLANG_RT.OSX_LIB} ; ${PETSC_TO_LIBRARY_LIB} ; ")
  set(PETSC_CLANG_RT.OSX_LIB "" CACHE STRING "" FORCE)
  set(PETSC_TO_LIBRARY_LIB "" CACHE STRING "" FORCE)
  message(STATUS "clang rt 2: ${PETSC_CLANG_RT.OSX_LIB} ; ${PETSC_TO_LIBRARY_LIB} ; ")
endif() 

FIND_PACKAGE( PETSc REQUIRED)
if ( NOT PETSC_FOUND )
  return()
endif()

#add_definitions( -DFEELPP_HAS_PETSC -DFEELPP_HAS_PETSC_H )
set(FEELPP_HAS_PETSC 1)
set(FEELPP_HAS_PETSC_H 1)

SET(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDES};${CMAKE_REQUIRED_INCLUDES}")
SET(FEELPP_LIBRARIES ${PETSC_LIBRARIES} ${FEELPP_LIBRARIES})

SET(BACKEND_PETSC petsc)
#INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR} ${PETSC_INCLUDE_CONF})
SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} PETSc/${PETSC_VERSION}" )

if (${PETSC_VERSION} VERSION_LESS 3.7)
  set(PETSC_CMAKE_CONFIG_PACKAGE_FILENAME PETScConfig.cmake)
else()
  set(PETSC_CMAKE_CONFIG_PACKAGE_FILENAME PETScBuildInternal.cmake)
endif()

find_path( PETSC_CMAKE_CONFIG_PACKAGE_DIR ${PETSC_CMAKE_CONFIG_PACKAGE_FILENAME}
  PATH
  ${PETSC_DIR}/conf
  ${PETSC_DIR}/lib/petsc/conf
  NO_DEFAULT_PATH )

#message("PETSC_CMAKE_CONG_PACKAGE_DIR=${PETSC_CMAKE_CONFIG_PACKAGE_DIR}")
if (NOT EXISTS ${PETSC_CMAKE_CONFIG_PACKAGE_DIR}/${PETSC_CMAKE_CONFIG_PACKAGE_FILENAME})
  macro (PETSC_TEST_CONF_MACRO petscmacro result)
    set (CMAKE_REQUIRED_INCLUDES ${PETSC_DIR}/include)
    # uncoment for force detection
    # if ( result )
    #  unset(result CACHE)
    #endif()
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <petscconf.h>
#ifdef ${petscmacro}
      int main() { return 0; }
#endif
 "
      ${result} )
    unset(CMAKE_REQUIRED_INCLUDES)
  endmacro()

  PETSC_TEST_CONF_MACRO( "PETSC_HAVE_MUMPS"  PETSC_HAVE_MUMPS )
  PETSC_TEST_CONF_MACRO( "PETSC_HAVE_PARMETISS" PETSC_HAVE_PARMETIS )
  PETSC_TEST_CONF_MACRO( "PETSC_HAVE_PTSCOTCH" PETSC_HAVE_PTSCOTCH )
  PETSC_TEST_CONF_MACRO( "PETSC_HAVE_ML" PETSC_HAVE_ML )
  PETSC_TEST_CONF_MACRO( "PETSC_HAVE_SUITESPARSE" PETSC_HAVE_SUITESPARSE )
  PETSC_TEST_CONF_MACRO( "PETSC_HAVE_HYPRE" PETSC_HAVE_HYPRE )
else()
  include(${PETSC_CMAKE_CONFIG_PACKAGE_DIR}/${PETSC_CMAKE_CONFIG_PACKAGE_FILENAME})

  list(APPEND FEELPP_LIBRARIES ${PETSC_PACKAGE_LIBS})
  #include_directories(${PETSC_PACKAGE_INCLUDES})
endif()

set(FEELPP_PETSC_ENABLED_OPTIONS)

if ( PETSC_GFORTRAN_LIB )
  set( GFORTRAN_LIBRARY ${PETSC_GFORTRAN_LIB} )
endif()

if ( PETSC_HAVE_MUMPS )
  set(FEELPP_HAS_MUMPS 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} MUMPS")
endif()

if ( PETSC_HAVE_PARMETIS )
  set(FEELPP_HAS_PARMETIS 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} ParMETIS")
endif()

if ( PETSC_HAVE_PTSCOTCH )
  set(FEELPP_HAS_SCOTCH 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} SCOTCH")
endif()

if ( PETSC_HAVE_ML )
  set(FEELPP_HAS_ML 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} ML")
endif()

if ( PETSC_HAVE_SUITESPARSE )
  set(FEELPP_HAS_SUITESPARSE 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} SUITESPARSE")
endif()

if ( PETSC_AMD_LIB )
  set(FEELPP_HAS_AMD_LIB 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} AMD")
endif()
if ( PETSC_COLAMD_LIB )
  set(FEELPP_HAS_COLAMD_LIB 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} COLAMD")
endif()
if ( PETSC_CHOLMOD_LIB )
  set(FEELPP_HAS_CHOLMOD_LIB 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} CHOLMOD")
endif()
if (PETSC_UMFPACK_LIB)
  set(FEELPP_HAS_UMFPACK_LIB 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} UMFPACK")
endif()

if ( PETSC_HAVE_HYPRE )
  set(FEELPP_HAS_HYPRE 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} HYPRE")
endif()

if ( PETSC_HAVE_YAML )
  set(FEELPP_HAS_YAML 1)
  set(FEELPP_PETSC_ENABLED_OPTIONS "${FEELPP_PETSC_ENABLED_OPTIONS} YAML")
endif()

if ( APPLE )
  message(STATUS "clang rt: ${PETSC_CLANG_RT.OSX_LIB} ; ${PETSC_TO_LIBRARY_LIB} ; ") 
  set(PETSC_CLANG_RT.OSX_LIB "" CACHE STRING "" FORCE)
  set(PETSC_TO_LIBRARY_LIB "" CACHE STRING "" FORCE)
  message(STATUS "clang rt 2: ${PETSC_CLANG_RT.OSX_LIB} ; ${PETSC_TO_LIBRARY_LIB} ; ") 
endif()
message( STATUS "Found PETSc with packages : ${FEELPP_PETSC_ENABLED_OPTIONS}")


  