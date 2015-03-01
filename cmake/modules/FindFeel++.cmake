# - Find Feel
# This module looks for Feel (Library for the Finite Element Method) support
# it will define the following values
#  FEELPP_INCLUDE_DIR = where feel/feelcore/feel.hpp can be found
#  FEELPP_LIBRARY    = the library to link in

# Check compiler
message(STATUS "[feelpp] Compiler version : ${CMAKE_CXX_COMPILER_ID}  ${CMAKE_CXX_COMPILER_VERSION}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # require at least gcc 4.7
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
    message(WARNING "GCC version must be at least 4.7!")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # require at least clang 3.3
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.3)
    message(WARNING "Clang version must be at least 3.3! we have clang ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    message(STATUS "[feelpp] Apple Clang version :  ${CMAKE_CXX_COMPILER_VERSION}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  message(STATUS "[feelpp] Intel version :  ${CMAKE_CXX_COMPILER_VERSION}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  message(STATUS "[feelpp] IBM XL compiler version :  ${CMAKE_CXX_COMPILER_VERSION}")
else()
  message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
  message(WARNING "CMAKE_CXX_COMPILER_ID=" ${CMAKE_CXX_COMPILER_ID})
endif()

#should check the version of gcc for -std=c++0x ou -std=c++11

IF( ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU") OR
    ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") OR
    ("${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang") OR
    ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel") )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -std=c++11" )
  if ( NOT ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel") )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=1024" )
  endif()
  if ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd3373" )
  endif()
  IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang")
    IF(APPLE OR FEELPP_USE_CLANG_LIBCXX)
      message(STATUS "[feelpp] Use clang libc++")
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++" )
    ELSE()
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++" )
    ENDIF()
  ENDIF()
ENDIF()

# IBM XL compiler
IF( ("${CMAKE_CXX_COMPILER_ID}" MATCHES "XL") )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qlanglvl=extc1x" )
endif()

LIST(REMOVE_DUPLICATES CMAKE_CXX_FLAGS)
LIST(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_DEBUG)
LIST(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_RELEASE)

INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckLibraryExists)

OPTION(FEELPP_ENABLE_SYSTEM_EIGEN3 "enable system eigen3 support" OFF)



OPTION(FEELPP_ENABLE_MOVE_SEMANTICS "enable move semantics(elision)" ON )
OPTION(FEELPP_ENABLE_INSTANTIATION_MODE "Instantiation mode" ON )
OPTION(FEELPP_ENABLE_MPI_MODE "Instantiation mode" ON )
OPTION(FEELPP_ENABLE_SCHED_SLURM "Enable Feel++ slurm submission scripts generation" OFF)
OPTION(FEELPP_ENABLE_SCHED_CCC "Enable Feel++ tgcc/ccc submission scripts generation" OFF)
OPTION(FEELPP_ENABLE_SCHED_LOADLEVELER "Enable Feel++ ibm(supermuc) submission scripts generation" OFF)
OPTION(FEELPP_ENABLE_TBB "enable feel++ TBB support" OFF)
OPTION(FEELPP_ENABLE_SLEPC "enable feel++ SLEPc support" ON)
OPTION(FEELPP_ENABLE_TRILINOS "enable feel++ Trilinos support" OFF)
OPTION(FEELPP_ENABLE_EXODUS "enable feel++ Exodus support" OFF)
if ( APPLE )
  OPTION(FEELPP_ENABLE_OPENTURNS "enable feel++ OpenTURNS support" OFF)
else()
  OPTION(FEELPP_ENABLE_OPENTURNS "enable feel++ OpenTURNS support" ON)
endif()
OPTION(FEELPP_ENABLE_OCTAVE "Enable Feel++/Octave interface" OFF)


OPTION(FEELPP_ENABLE_OPENGL "enable feel++ OpenGL support" ON)
OPTION(FEELPP_DISABLE_EIGEN_ALIGNMENT "disable alignement (hence vectorization) in Eigen" OFF)

# enable mpi mode
IF ( FEELPP_ENABLE_MPI_MODE )
  SET( FEELPP_ENABLE_MPI_MODE 1 )
ENDIF()

# disable alignement
MARK_AS_ADVANCED(FEELPP_DISABLE_EIGEN_ALIGNMENT)
if ( FEELPP_DISABLE_EIGEN_ALIGNMENT )
  add_definitions(-DEIGEN_DONT_ALIGN=1 -DEIGEN_DONT_VECTORIZE=1)
  message(STATUS "[feelpp] Disabling alignment and vectorisation in Feel++/Eigen")
endif()

# enable move semantics
MARK_AS_ADVANCED(FEELPP_ENABLE_MOVE_SEMANTICS)

# enable instantiation
MARK_AS_ADVANCED(FEELPP_ENABLE_INSTANTIATION_MODE)
IF ( FEELPP_ENABLE_INSTANTIATION_MODE )
  SET( FEELPP_INSTANTIATION_MODE 1 )
ENDIF()
SET(FEELPP_MESH_MAX_ORDER "5" CACHE STRING "maximum geometrical order in templates to instantiate" )

# enable host specific
include(feelpp.host)

if ( FEELPP_ENABLE_TBB )
  FIND_PACKAGE(TBB)
  IF ( TBB_FOUND )
    INCLUDE_DIRECTORIES( ${TBB_INCLUDE_DIR} )
    SET(FEELPP_LIBRARIES ${TBB_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Tbb" )
  ENDIF (TBB_FOUND )
endif()

# only activate OpenMP for gcc
# (clang support should be ok by 3.5)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    OPTION( FEELPP_ENABLE_OPENMP "Enable OpenMP" OFF )
    if ( FEELPP_ENABLE_OPENMP )
        FIND_PACKAGE(OpenMP)

        if(OPENMP_FOUND)
            set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
            set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
            SET( FEELPP_HAS_OPENMP 1 )
            SET( FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} OpenMP" )
        endif()
    endif()
endif()

# on APPLE enfore the use of macports openmpi version
if ( APPLE )
  if ( EXISTS /opt/local/lib/openmpi/bin/mpic++ )
    set(MPI_COMPILER /opt/local/lib/openmpi/bin/mpic++)
  endif()

  #  set(MPI_LIBRARY "MPI_LIBRARY-NOTFOUND" )
  MESSAGE(STATUS "[feelpp] Use mpi compiler ${MPI_COMPILER}")

endif( APPLE )
FIND_PACKAGE(MPI REQUIRED)
IF ( MPI_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")
  SET( FEELPP_HAS_MPI 1 )
  SET( FEELPP_HAS_MPI_H 1 )
  ADD_DEFINITIONS( -DFEELPP_HAS_MPI=1 -DFEELPP_HAS_MPI_H=1 )
  SET(FEELPP_BOOST_MPI mpi)
  SET(FEELPP_LIBRARIES ${MPI_LIBRARIES} ${FEELPP_LIBRARIES})
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Mpi" )

  # Check for MPI IO Support

  #TRY_COMPILE(MPIIO_SUCCESS ${CMAKE_CURRENT_BINARY_DIR}/tryCompileMPIIO
  #${CMAKE_SOURCE_DIR}/cmake/codes/try-mpiio.cpp
  #LINK_LIBRARIES ${FEELPP_LIBRARIES} )
  set(CMAKE_REQUIRED_LIBRARIES_save ${CMAKE_REQUIRED_LIBRARIES})
  set(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})
  set(CMAKE_REQUIRED_INCLUDES_save ${CMAKE_REQUIRED_INCLUDES})
  set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
  CHECK_CXX_SOURCE_COMPILES(
      "
      #include <mpi.h>

      int main(int argc, char** argv)
      {
      MPI_File fh;
      MPI_Status status;
      MPI_Info info;
      }
      "
      MPIIO_DETECTED)

  # Check if we have the types from the 2.2 standard
  # needed for MPI IO
  CHECK_CXX_SOURCE_COMPILES(
      "
      #include <mpi.h>

      int main(int argc, char** argv)
      {
      MPI_INT32_T i32;
      MPI_INT64_T i64;
      }
      "
      MPIIO_HAS_STD_22_TYPES)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_save})
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES_save})

  # should be compatible with 2.2 standard
  IF ( MPIIO_DETECTED AND NOT MPIIO_HAS_STD_22_TYPES)
      include(CheckTypeSize)
      check_type_size("int" SIZEOF_INT BUILTIN_TYPES_ONLY)
      IF(SIZEOF_INT STREQUAL 4)
          SET(FEELPP_MPI_INT32 MPI_INT)
      ELSE()
          MESSAGE(STATUS "[feelpp] MPIIO: Cannot find a compatible int32 type")
          MESSAGE(STATUS "[feelpp] MPIIO: disabling MPIIO")
          set(MPIIO_DETECTED 0)
      ENDIF()
      check_type_size("long" SIZEOF_LONG_LONG BUILTIN_TYPES_ONLY)
      IF(SIZEOF_LONG_LONG STREQUAL 8)
          SET(FEELPP_MPI_INT64 MPI_LONG_LONG)
      ELSE()
         MESSAGE(STATUS "[feelpp] MPIIO: Cannot find a compatible int64 type")
         MESSAGE(STATUS "[feelpp] MPIIO: disabling MPIIO")
         set(MPIIO_DETECTED 0)
      ENDIF()
  ENDIF()

  IF (MPIIO_DETECTED)
      MESSAGE(STATUS "[feelpp] MPIIO detected and enabled.")
      SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Mpi-IO" )
      SET(FEELPP_HAS_MPIIO 1)
  ELSE()
      MESSAGE(WARNING "MPIIO not detected and disabled (Related features disable, e.g. Ensight Gold exporter).")
  ENDIF()
ENDIF()

CHECK_FUNCTION_EXISTS(fmemopen FEELPP_HAS_STDIO_FMEMOPEN)
MESSAGE(STATUS "[feelpp] FMemOpen: ${FEELPP_HAS_STDIO_FMEMOPEN}")


Check_Include_File_CXX(dlfcn.h FEELPP_HAS_DLFCN_H)
if ( FEELPP_HAS_DLFCN_H )
  add_definitions(-DFEELPP_HAS_DLFCN_H)
endif()
CHECK_LIBRARY_EXISTS (dl dlopen "" FEELPP_HAS_LIBDL)
IF (FEELPP_HAS_LIBDL)
  set(DL_LIBRARY dl)
  SET (FEELPP_HAS_DLOPEN 1)
  add_definitions(-DFEELPP_HAS_DLOPEN)
  SET(FEELPP_LIBRARIES ${DL_LIBRARY} ${FEELPP_LIBRARIES})
ELSE ()
  CHECK_FUNCTION_EXISTS (dlopen FEELPP_HAS_DLOPEN)
ENDIF (FEELPP_HAS_LIBDL)

find_package(GMP)
if ( GMP_FOUND )
  SET(FEELPP_LIBRARIES  ${GMP_LIBRARIES} ${FEELPP_LIBRARIES})
  message(STATUS "[feelpp] GMP: ${GMP_LIBRARIES}" )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Gmp" )
endif()

#
# Intel MKL
# 
OPTION( FEELPP_ENABLE_MKL "Enable the Intel MKL library" OFF )
if ( FEELPP_ENABLE_MKL )
  find_package(MKL)
  if ( MKL_FOUND )
    
    message(STATUS "[feelpp] MKL Includes: ${MKL_INCLUDE_DIRS}")
    message(STATUS "[feelpp] MKL Libraries: ${MKL_LIBRARIES}")
    set(FEELPP_HAS_MKL 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Intel(MKL)" )
    include_directories( ${MKL_INCLUDE_DIRS} )
    SET(FEELPP_LIBRARIES ${MKL_LIBRARIES} ${FEELPP_LIBRARIES})
    #  enable MKL wherever possible for eigen3
    add_definitions(-DEIGEN_USE_MKL_ALL=1)
  else( MKL_FOUND )
    #
    # Blas and Lapack
    #
    if (APPLE)
      # FIND_LIBRARY(ATLAS_LIBRARY
      #   NAMES
      #   atlas
      #   PATHS
      #   /opt/local/lib/lib
      #   NO_DEFAULT_PATH
      #   )
      # message(STATUS "[feelpp] ATLAS: ${ATLAS_LIBRARY}" )
      # IF( ATLAS_LIBRARY )
      #   SET(FEELPP_LIBRARIES ${ATLAS_LIBRARY} ${FEELPP_LIBRARIES})
      # ENDIF()
      FIND_PACKAGE(LAPACK)
    else (APPLE)
      FIND_PACKAGE(LAPACK)
    endif (APPLE)
    SET(FEELPP_LIBRARIES  ${LAPACK_LIBRARIES} ${FEELPP_LIBRARIES})
  endif(MKL_FOUND)
endif(FEELPP_ENABLE_MKL)

# HDF5
# On debian, 
# - do not install hdf5-helpers, otherwise it will pick the serial version by default
# - install only the libhdf5-openmpi-dev package
OPTION( FEELPP_ENABLE_HDF5 "Enable the HDF5 library" ON )
if ( FEELPP_ENABLE_HDF5 )
  find_package(HDF5)
  if( HDF5_FOUND ) 
    if( HDF5_IS_PARALLEL )
        message(STATUS "[feelpp] HDF5 - Headers ${HDF5_INCLUDE_DIRS}" )
        message(STATUS "[feelpp] HDF5 - Libraries ${HDF5_LIBRARIES}" )
        include_directories( ${HDF5_INCLUDE_DIRS} )
        set(FEELPP_LIBRARIES ${HDF5_LIBRARIES} ${FEELPP_LIBRARIES})
        set(FEELPP_HAS_HDF5 1)
        set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} HDF5" )
    else()
        MESSAGE(STATUS "[feelpp] HDF5 has been found but is not parallel, HDF5 is not enabled in Feel++")
    endif()
  endif()
endif(FEELPP_ENABLE_HDF5)


# XDMF
find_package(XDMF QUIET)
if (XDMF_FOUND)
    include_directories( ${XDMF_INCLUDE_DIRS} )
    set(FEELPP_LIBRARIES ${XDMF_LIBRARIES})
    set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} XDMF" )
    message(STATUS "Found Xdmf." )
else()
    message(STATUS "Could not find Xdmf." )
endif (XDMF_FOUND)

option(FEELPP_ENABLE_PYTHON_WRAPPING "Enable Boost.Python wrapping implementation" OFF)

# Boost
SET(BOOST_MIN_VERSION "1.55.0")

# Making consecutive calls to find_package for Boost to find optional components (boost_python for now)
# Making only one call to find_package and having one of the component not installed will mark Boost as not found

# First we try to find boost with the python components
FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} COMPONENTS python )
if(Boost_PYTHON_FOUND)
    set(FEELPP_HAS_BOOST_PYTHON 1)
    set(FEELPP_LIBRARIES ${Boost_PYTHON_LIBRARY} ${FEELPP_LIBRARIES})

    if(FEELPP_ENABLE_PYTHON_WRAPPING)
        SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Python-Wrapping" )
    endif()
else()
    if(FEELPP_ENABLE_PYTHON_WRAPPING)
        message(FATAL_ERROR "[feelpp] Boost.Python was not found on your system (Required for Python Wrapping)." )
    else()
        message(STATUS "[feelpp] Boost.Python was not found on your system." )
    endif()
endif()

# Then we try to find rest of the Boost components
FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} REQUIRED date_time filesystem system program_options unit_test_framework signals ${FEELPP_BOOST_MPI} regex serialization )
if(Boost_FOUND)
  IF(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "51")
    add_definitions(-DBOOST_RESULT_OF_USE_TR1)
    message(STATUS "[feelpp] added -DBOOST_RESULT_OF_USE_TR1" )
  ENDIF()
  IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    # ensures that boost.signals2 compiles with clang++ >= 3.1
    IF(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "52")
      #add_definitions(-DBOOST_NO_CXX11_VARIADIC_TEMPLATES)
      #message(STATUS "[feelpp] added -DBOOST_NO_CXX11_VARIADIC_TEMPLATES" )
    ELSE()
      add_definitions(-DBOOST_NO_VARIADIC_TEMPLATES)
      message(STATUS "[feelpp] added -DBOOST_NO_VARIADIC_TEMPLATES" )
    ENDIF()

    # This is meant to solves a compilation issue arising with Boost 1.56 and clang compilers (confirmed for 3.4 and 3.5)
    IF(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "55")
      add_definitions(-DBOOST_PP_VARIADICS=0)
      message(STATUS "[feelpp] added -DBOOST_PP_VARIADICS=0" )
    ENDIF()
  ENDIF()
else()
  message(STATUS "[feelpp] Please check your boost version - Should be at least ${BOOST_MIN_VERSION}")
endif()

IF ( FEELPP_ENABLE_MOVE_SEMANTICS AND Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION LESS "57" )
  SET( BOOST_UBLAS_MOVE_SEMANTICS 1 CACHE STRING "Enable Boost Ublas move semantics" FORCE )
  ADD_DEFINITIONS( -DBOOST_UBLAS_MOVE_SEMANTICS )
ENDIF()


OPTION(BOOST_ENABLE_TEST_DYN_LINK "enable boost test with dynamic lib" ON)
MARK_AS_ADVANCED(BOOST_ENABLE_TEST_DYN_LINK)

set(Boost_ADDITIONAL_VERSIONS "1.39" "1.40" "1.41" "1.42" "1.43" "1.44" "1.45" "1.46" "1.47" "1.48" "1.49" "1.50" "1.51" "1.52" "1.53" "1.54" "1.55")
set( BOOST_PARAMETER_MAX_ARITY 24 )
#set( BOOST_FILESYSTEM_VERSION 2)
set( BOOST_FILESYSTEM_VERSION 3)
if (BOOST_ENABLE_TEST_DYN_LINK)
  add_definitions( -DBOOST_PARAMETER_MAX_ARITY=${BOOST_PARAMETER_MAX_ARITY} -DBOOST_TEST_DYN_LINK -DBOOST_FILESYSTEM_VERSION=${BOOST_FILESYSTEM_VERSION})
else (BOOST_ENABLE_TEST_DYN_LINK)
  add_definitions( -DBOOST_PARAMETER_MAX_ARITY=${BOOST_PARAMETER_MAX_ARITY} -DBOOST_FILESYSTEM_VERSION=${BOOST_FILESYSTEM_VERSION})
endif (BOOST_ENABLE_TEST_DYN_LINK)

# undefined BOOST_UBLAS_TYPE_CHECK
add_definitions(-UBOOST_UBLAS_TYPE_CHECK )

# this fix an issue with boost filesystem: boost is usually no compiled with
# std=c++0x and we compile with it, this causes problems with the macro
# BOOST_SCOPED_ENUM macros whose behavior differs in both case and would
# generate different c++ codes and undefined references at link time.
# in a short future, this should not be necessary anymore
IF(NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR NOT APPLE)
  ADD_DEFINITIONS(-DBOOST_NO_SCOPED_ENUMS)
  IF(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "51")
    ADD_DEFINITIONS(-DBOOST_NO_CXX11_SCOPED_ENUMS)
  endif()
endif()

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR}   ${BOOST_INCLUDE_PATH})

SET(FEELPP_LIBRARIES ${Boost_LIBRARIES} ${FEELPP_LIBRARIES})

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/feel)

INCLUDE_DIRECTORIES(BEFORE contrib/)


#FIND_PACKAGE(GINAC)
#IF( GINAC_FOUND )
#  set( FEELPP_HAS_GINAC 1 )
#  INCLUDE_DIRECTORIES( GINAC_INCLUDE_DIRS )
#  SET(FEELPP_LIBRARIES ${GINAC_LIBRARIES} ${FEELPP_LIBRARIES})
#  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GINAC" )
#ENDIF()

add_definitions(-DHAVE_LIBDL)

OPTION( FEELPP_ENABLE_NT2 "Enable the numerical toolkit tmplate library" OFF )
if ( FEELPP_ENABLE_NT2 )
  #set(NT2_SOURCE_ROOT ${FEELPP_ROOT}/contrib/nt2)
  #set(NT2_WITH_TESTS OFF)
  option(NT2_WITH_TESTS "Enable benchmarks and unit tests" OFF)
  add_subdirectory(contrib/nt2)
  set(CMAKE_CXX_FLAGS "${NT2_SIMD_FLAGS} ${CMAKE_CXX_FLAGS}")
  include_directories(${CMAKE_BINARY_DIR}/contrib/nt2/include/)

  foreach(module ${NT2_FOUND_COMPONENTS})
    string(TOUPPER ${module} module_U)
    if(NT2_${module_U}_ROOT)
      include_directories(${NT2_${module_U}_ROOT}/include)
      message(status "[feelpp/nt2] adding ${NT2_${module_U}_ROOT}/include" )
    endif()
  endforeach()

  set(NT2_FOUND 1)
  set(FEELPP_HAS_NT2 1)
  SET(FEELPP_LIBRARIES nt2  ${FEELPP_LIBRARIES}  )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} NT2" )
  message(STATUS "[feelpp] nt2 is enabled" )
  
endif( FEELPP_ENABLE_NT2 )

if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/feel AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/contrib )

  #
  # cln and ginac
  #
  add_definitions(-DIN_GINAC -DHAVE_LIBDL)
  include_directories(${FEELPP_BUILD_DIR}/contrib/cln/include ${FEELPP_SOURCE_DIR}/contrib/ginac/ ${FEELPP_BUILD_DIR}/contrib/ginac/ ${FEELPP_SOURCE_DIR}/contrib/ginac/ginac ${FEELPP_BUILD_DIR}/contrib/ginac/ginac )
  SET(FEELPP_LIBRARIES feelpp_ginac ${CLN_LIBRARIES} ${FEELPP_LIBRARIES} ${CMAKE_DL_LIBS} )
  set(DL_LIBS ${CMAKE_DL_LIBS})
  add_subdirectory(contrib/ginac)

endif()

#
# submodules
#
include(feelpp.module.hpddm)
include(feelpp.module.nlopt)
include(feelpp.module.ipopt)
include(feelpp.module.cereal)

#
# HARTS
#
OPTION( FEELPP_ENABLE_HARTS "Enable Harts (Runtime parallelization system)" OFF )
if ( FEELPP_ENABLE_HARTS )
  FIND_PACKAGE( HARTS )
  if( HARTS_FOUND )
    SET(CMAKE_REQUIRED_INCLUDES ${HARTS_INCLUDES} ${CMAKE_REQUIRED_INCLUDES})
    INCLUDE_DIRECTORIES( ${HARTS_INCLUDES} )
    ADD_DEFINITIONS(${HARTS_DEFINITIONS})
    SET(FEELPP_LIBRARIES ${HARTS_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(FEELPP_HAS_HARTS 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} HARTS" )
  endif()
endif()


if ( FEELPP_ENABLE_EXODUS )
  include_directories(${FEELPP_SOURCE_DIR}/contrib/exodus-5.24/exodus/cbind/include/)
  add_subdirectory(contrib/exodus-5.24/exodus)
  #add_subdirectory(contrib/exodus-5.24/nemesis)
  set(FEELPP_HAS_EXODUS 1)
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Exodus" )
endif()

#
# Eigen
#
if ( FEELPP_ENABLE_SYSTEM_EIGEN3 )
  FIND_PACKAGE(Eigen3)
endif()
if (NOT EIGEN3_FOUND AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/feel AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/contrib )
  option(EIGEN_BUILD_PKGCONFIG "Build pkg-config .pc file for Eigen" OFF)
  set(EIGEN_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/feel)
  add_subdirectory(contrib/eigen)
  set( EIGEN3_INCLUDE_DIR ${FEELPP_SOURCE_DIR}/contrib/eigen )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Eigen3/Contrib" )
elseif( EIGEN3_FOUND )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Eigen3/System" )
endif()
INCLUDE_DIRECTORIES( ${EIGEN3_INCLUDE_DIR} )
message(STATUS "[feelpp] Eigen3: ${EIGEN3_INCLUDE_DIR}" )

#FIND_PACKAGE(Eigen2 REQUIRED)
#INCLUDE_DIRECTORIES( ${Eigen2_INCLUDE_DIR} )
#add_subdirectory(contrib/eigen)
#INCLUDE_DIRECTORIES( ${FEELPP_SOURCE_DIR}/contrib/eigen )
#add_definitions( -DEIGEN_NO_STATIC_ASSERT )

#
# Metis
#
FIND_PACKAGE(Metis)
if ( METIS_FOUND )
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
  #  LINK_DIRECTORIES(${METIS_LIBRARIES})
  SET(FEELPP_LIBRARIES ${METIS_LIBRARY} ${FEELPP_LIBRARIES})
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Metis" )
endif( METIS_FOUND )

#
# Ann
#
FIND_PACKAGE(ANN)
if ( ANN_FOUND )
  set(FEELPP_HAS_ANN_H 1)
  INCLUDE_DIRECTORIES( ${ANN_INCLUDE_DIR} )
  SET(FEELPP_LIBRARIES ${ANN_LIBRARIES} ${FEELPP_LIBRARIES})
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} ANN" )
endif()

#
# GLPK
#
FIND_PACKAGE(GLPK)
if ( GLPK_FOUND )
  set(FEELPP_HAS_GLPK_H 1)
  INCLUDE_DIRECTORIES( ${GLPK_INCLUDE_DIR} )
  SET(FEELPP_LIBRARIES ${GLPK_LIBRARIES} ${FEELPP_LIBRARIES})
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GLPK" )
endif()

# google perf tools
option(FEELPP_ENABLE_GOOGLEPERFTOOLS "Enable Google Perf Tools (tcmalloc, stracktrace and profiler)" OFF)
if ( FEELPP_ENABLE_GOOGLEPERFTOOLS )
  find_package(GooglePerfTools)
  if ( GOOGLE_PERFTOOLS_FOUND )
    set(FEELPP_HAS_GPERFTOOLS 1 )
    message(STATUS "[feelpp] Google PerfTools: ${TCMALLOC_LIBRARIES} ${STACKTRACE_LIBRARIES} ${PROFILER_LIBRARIES}")
    include_directories(${GOOGLE_PERFTOOLS_INCLUDE_DIR})
    SET(FEELPP_LIBRARIES  ${FEELPP_LIBRARIES} ${TCMALLOC_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GooglePerfTools" )
  endif()
endif( FEELPP_ENABLE_GOOGLEPERFTOOLS )

option(FEELPP_ENABLE_DDT "Enable DDT support" OFF)
if ( FEELPP_ENABLE_DDT )
  find_package(DDT)
  if ( DDT_FOUND )
    message(STATUS "[feelpp] DDT: ${DDT_LIBRARIES}")
    SET(FEELPP_LIBRARIES  ${FEELPP_LIBRARIES} ${DDT_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} DDT" )
  endif()
endif( FEELPP_ENABLE_DDT )

# google gflags
find_package(GFLAGS REQUIRED)

INCLUDE_DIRECTORIES( ${GFLAGS_INCLUDE_DIR} )

SET(FEELPP_LIBRARIES ${GFLAGS_LIBRARIES} ${FEELPP_LIBRARIES})
if ( ${GFLAGS_INCLUDE_DIR} MATCHES "/contrib/" )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GFLAGS/Contrib" )
else()
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GFLAGS/System" )
endif()

# google glog
find_package(GLOG REQUIRED)

INCLUDE_DIRECTORIES( ${GLOG_INCLUDE_DIR} )
SET(FEELPP_LIBRARIES ${GLOG_LIBRARIES} ${FEELPP_LIBRARIES})
if ( ${GLOG_INCLUDE_DIR} MATCHES "/contrib/" )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GLOG/Contrib" )
else()
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GLOG/System" )
endif()


# xml
find_package(LibXml2 2.6.27)
if ( LIBXML2_FOUND )
    message(STATUS "[feelpp] LibXml2: ${LIBXML2_INCLUDE_DIR} ${LIBXML2_LIBRARIES}")
    include_directories(${LIBXML2_INCLUDE_DIR})
    SET(FEELPP_LIBRARIES ${LIBXML2_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} LibXml2" )
    set( FEELPP_HAS_LIBXML2 1 )
endif()

# Python libs
FIND_PACKAGE(PythonLibs)
if ( PYTHONLIBS_FOUND )
  message(STATUS "[feelpp] PythonLibs: ${PYTHON_INCLUDE_DIRS} ${PYTHON_LIBRARIES}")
  include_directories(${PYTHON_INCLUDE_DIRS})
  SET(FEELPP_LIBRARIES ${PYTHON_LIBRARIES} ${FEELPP_LIBRARIES})
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Python" )
  set( FEELPP_HAS_PYTHON 1 )
endif()

#
# Python interp
#
FIND_PACKAGE(PythonInterp REQUIRED)
if(PYTHONINTERP_FOUND)
  execute_process(COMMAND
    ${PYTHON_EXECUTABLE}
    -c "import sys; print sys.version[0:3]"
    OUTPUT_VARIABLE PYTHON_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  message(STATUS "[feelpp] Found python version ${PYTHON_VERSION}")
endif()

# metis
FIND_LIBRARY(METIS_LIBRARY
  NAMES
  metis
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  #    "/opt/local/lib"
  )
message(STATUS "[feelpp] Metis: ${METIS_LIBRARY}" )
IF( METIS_LIBRARY )
  SET(FEELPP_LIBRARIES ${METIS_LIBRARY} ${FEELPP_LIBRARIES})
ENDIF()

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
ENDIF()

FIND_PACKAGE(Scotch)
IF( SCOTCH_FOUND )
  message(STATUS "[feelpp] SCOTCH: ${SCOTCH_LIBRARIES}" )
  SET(FEELPP_LIBRARIES ${SCOTCH_LIBRARIES} ${FEELPP_LIBRARIES})
ENDIF()

find_package(ML)
message(STATUS "[feelpp] ML: ${ML_LIBRARY}" )
IF ( ML_FOUND )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} PETSc" )
  SET(FEELPP_LIBRARIES ${ML_LIBRARY} ${FEELPP_LIBRARIES})
ENDIF()

if ( NOT GFORTRAN_LIBRARY )
  FIND_LIBRARY(GFORTRAN_LIBRARY
    NAMES
    gfortran
    PATHS
    $ENV{LIBRARY_PATH}
    /opt/local/lib
    /usr/lib/gcc/x86_64-linux-gnu/
    PATH_SUFFIXES
    gcc47 gcc46 gcc45 gcc44 4.7 4.6 4.5 4.4
    )
endif()

message(STATUS "[feelpp] gfortran lib: ${GFORTRAN_LIBRARY} ")
if ( GFORTRAN_LIBRARY )
  set( FEELPP_LIBRARIES ${GFORTRAN_LIBRARY} ${FEELPP_LIBRARIES})
endif()
FIND_PACKAGE(MUMPS)
if ( GFORTRAN_LIBRARY AND MUMPS_FOUND )
  set( FEELPP_HAS_MUMPS 1 )
  set( FEELPP_LIBRARIES ${MUMPS_LIBRARIES} ${FEELPP_LIBRARIES} )
endif()
FIND_LIBRARY(SUITESPARSECONFIG_LIBRARY
  NAMES
  suitesparseconfig
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  $ENV{SUITESPARSE_DIR}/lib
  )
IF ( SUITESPARSECONFIG_LIBRARY )
  SET(FEELPP_LIBRARIES  ${SUITESPARSECONFIG_LIBRARY} ${FEELPP_LIBRARIES})
endif()
FIND_LIBRARY(AMD_LIBRARY
  NAMES
  amd
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  $ENV{SUITESPARSE_DIR}/lib
  )
IF ( AMD_LIBRARY )
  SET(FEELPP_LIBRARIES  ${AMD_LIBRARY} ${FEELPP_LIBRARIES})
endif()

FIND_LIBRARY(COLAMD_LIBRARY
  NAMES
  colamd
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  $ENV{SUITESPARSE_DIR}/lib
  )
IF ( COLAMD_LIBRARY )
  SET(FEELPP_LIBRARIES  ${COLAMD_LIBRARY} ${FEELPP_LIBRARIES})
endif()

FIND_LIBRARY(CHOLMOD_LIBRARY
  NAMES
  cholmod
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  $ENV{SUITESPARSE_DIR}/lib
  )

FIND_LIBRARY(UMFPACK_LIBRARY
  NAMES
  umfpack
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  $ENV{SUITESPARSE_DIR}/lib
  )
message(STATUS "[feelpp] SuiteSparseConfig: ${SUITESPARSECONFIG_LIBRARY}" )
message(STATUS "[feelpp] Amd: ${AMD_LIBRARY}" )
message(STATUS "[feelpp] ColAmd: ${COLAMD_LIBRARY}" )
message(STATUS "[feelpp] Cholmod: ${CHOLMOD_LIBRARY}" )
message(STATUS "[feelpp] Umfpack: ${UMFPACK_LIBRARY}" )
if ( AMD_LIBRARY AND CHOLMOD_LIBRARY AND UMFPACK_LIBRARY )
  SET(FEELPP_LIBRARIES ${UMFPACK_LIBRARY} ${CHOLMOD_LIBRARY} ${FEELPP_LIBRARIES})
endif()

FIND_LIBRARY(YAML_LIBRARY
  NAMES
  yaml
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  $ENV{SUITESPARSE_DIR}/lib
  /opt/local/lib
  )
if ( YAML_LIBRARY )
  SET(FEELPP_LIBRARIES ${YAML_LIBRARY} ${FEELPP_LIBRARIES})
endif()

#
# Petsc
#
FIND_PACKAGE( PETSc REQUIRED)
if ( PETSC_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDES};${CMAKE_REQUIRED_INCLUDES}")
  SET(FEELPP_LIBRARIES ${PETSC_LIBRARIES} ${FEELPP_LIBRARIES})
  SET(BACKEND_PETSC petsc)
  INCLUDE_DIRECTORIES(
    ${PETSC_INCLUDE_DIR}
    ${PETSC_INCLUDE_CONF}
    )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} PETSc" )

endif( PETSC_FOUND )

# ML was already searched for, if it was not found then try again to look for it
# in PETSC_DIR
if ( NOT ML_FOUND )
  find_package(ML)
  message(STATUS "[feelpp] ML(PETSc): ${ML_LIBRARY}" )
  IF ( ML_LIBRARY )
    SET(FEELPP_LIBRARIES ${ML_LIBRARY} ${FEELPP_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} ML" )
  ENDIF()
endif()

#
# parpack
#
FIND_LIBRARY(PARPACK_LIBRARY NAMES parpack)
if (PARPACK_LIBRARY)
  SET(PARPACK_LIBRARIES ${PARPACK_LIBRARY})
  SET(FEELPP_LIBRARIES ${PARPACK_LIBRARIES} ${FEELPP_LIBRARIES})
endif()
MARK_AS_ADVANCED( PARPACK_LIBRARY )


#
# SLEPc
#
if (FEELPP_ENABLE_SLEPC)
  FIND_PACKAGE( SLEPc )
  if ( SLEPC_FOUND )
    SET(CMAKE_REQUIRED_INCLUDES "${SLEPC_INCLUDES};${CMAKE_REQUIRED_INCLUDES}")
    INCLUDE_DIRECTORIES( ${SLEPC_INCLUDE_DIR} )
    SET(FEELPP_LIBRARIES ${SLEPC_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} SLEPc" )
  endif(SLEPC_FOUND)
endif(FEELPP_ENABLE_SLEPC)


#
# Trilinos
#
if (FEELPP_ENABLE_TRILINOS)
  FIND_PACKAGE(Trilinos)
  if ( TRILINOS_FOUND )
    INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIR})
    SET(FEELPP_LIBRARIES ${TRILINOS_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(BACKEND_TRILINOS trilinos)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Trilinos" )
  endif( TRILINOS_FOUND )
endif (FEELPP_ENABLE_TRILINOS)

#
# OpenTURNS
#
IF ( FEELPP_ENABLE_OPENTURNS )
  FIND_PACKAGE( OpenTURNS )
  if ( OPENTURNS_FOUND )
    MESSAGE(STATUS "[feelpp] OpenTURNS Libraries: ${OpenTURNS_LIBRARIES}")
    MESSAGE(STATUS "[feelpp] OpenTURNS Headers: ${OpenTURNS_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES(${OpenTURNS_INCLUDE_DIRS})
    #SET(FEELPP_LIBRARIES ${OpenTURNS_LIBRARIES} ${FEELPP_LIBRARIES})
    # now OpenTURNS_LIBRARIES are used in crb_add_python_module
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} OpenTURNS" )
  endif( OPENTURNS_FOUND )
endif()

#
# VTK
#
OPTION( FEELPP_ENABLE_VTK "Enable the VTK library" ON )
OPTION( FEELPP_ENABLE_VTK_INSITU "Enable In-Situ Visualization using VTK/Paraview" OFF )
if ( FEELPP_ENABLE_VTK )

    # If we enable in-situ visualization
    # We need to look for the Paraview package for the corresponding headers
    # As Paravie integrates vtk headers we don't need them
    if ( FEELPP_ENABLE_VTK_INSITU )
        FIND_PACKAGE(ParaView REQUIRED COMPONENTS vtkParallelMPI vtkPVCatalyst vtkPVPythonCatalyst)
        message(STATUS "[ParaView] Use file: ${PARAVIEW_USE_FILE}")
        INCLUDE(${PARAVIEW_USE_FILE})

        # Enable VTK exporter and insitu in config
        set(FEELPP_VTK_INSITU_ENABLED 1)

        # Mark VTK as available
        set(FEELPP_HAS_VTK 1)
        # Check for version to ensure that we are able to
        # use an external communicator
        set(VTK_HAS_PARALLEL 0)
        if( VTK_MAJOR_VERSION EQUAL 6 OR VTK_MAJOR_VERSION GREATER 6 )
            set(VTK_HAS_PARALLEL 1)
        endif()

        SET(FEELPP_LIBRARIES ${ParaView_LIBRARIES} ${FEELPP_LIBRARIES})
        SET(FEELPP_LIBRARIES ${VTK_LIBRARIES} ${FEELPP_LIBRARIES})
        SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} ParaView/VTK" )

        message(STATUS "Found ParaView ${PARAVIEW_VERSION_FULL}/VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
    else()
        FIND_PACKAGE(VTK)
        if( VTK_FOUND )
            set(FEELPP_HAS_VTK 1)
            MESSAGE(STATUS "[feelpp] Found VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")# ${VTK_LIBRARIES}")

            # Check for MPI suppot in VTK
            set(VTK_HAS_PARALLEL 0)
            # Prior to VTK version 6, VTK_KITS was used
            if( VTK_MAJOR_VERSION LESS 6 )
                #message("Available VTK KITS: ${VTK_KITS}")
                list(FIND VTK_KITS "PARALLEL" __test_vtk_parallel)
                if( NOT ( ${__test_vtk_parallel} EQUAL -1 ) )
                    set(VTK_HAS_PARALLEL 1)
                    message(WARNING "External initialization of MPI Communicator is not available in VTK 5. VTK parallel export will be disabled.")
                else()
                    message(WARNING "MPI support for VTK 5 is not activated. VTK parallel export will be disabled.")
                endif()
                unset(__test_vtk_parallel)
                # From version 6 modules ared used
            else()
                #message(FATAL_ERROR "${VTK_MODULES_ENABLED}")
                list(FIND VTK_MODULES_ENABLED "vtkParallelMPI" __test_vtk_parallel)
                if( NOT (${__test_vtk_parallel} EQUAL -1) )
                    set(VTK_HAS_PARALLEL 1)
                else()
                    message(WARNING "MPI support for VTK 6 is not activated. VTK parallel export will be disabled.")
                endif()
                unset(__test_vtk_parallel)
            endif() 

            if ( NOT FEELPP_ENABLE_OPENGL )
                SET(VTK_LIBRARIES "-lvtkRendering -lvtkGraphics -lvtkImaging  -lvtkFiltering -lvtkCommon -lvtksys" )
            endif()
            INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
            MARK_AS_ADVANCED( VTK_DIR )
            SET(FEELPP_LIBRARIES ${VTK_LIBRARIES} ${FEELPP_LIBRARIES})
            SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} VTK" )

        endif()
    endif()
endif( FEELPP_ENABLE_VTK )

#
# Octave
#
if ( FEELPP_ENABLE_OCTAVE )
  FIND_PACKAGE(Octave)
  if ( OCTAVE_FOUND )

    # find octave-config and get variables from it
    FIND_PROGRAM(OCTAVE_CONFIG octave-config)
    IF(NOT OCTAVE_CONFIG)
      MESSAGE(FATAL_ERROR "Octave is required, but octave-config was not found.  Please install Octave and rerun cmake.")
    ENDIF(NOT OCTAVE_CONFIG)

    EXECUTE_PROCESS(COMMAND ${OCTAVE_CONFIG} --oct-site-dir
      OUTPUT_VARIABLE OCT_SITE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    EXECUTE_PROCESS(COMMAND ${OCTAVE_CONFIG} --m-site-dir
      OUTPUT_VARIABLE M_SITE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    EXECUTE_PROCESS(COMMAND ${OCTAVE_CONFIG} -p OCTINCLUDEDIR
      OUTPUT_VARIABLE OCTINCLUDEDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    EXECUTE_PROCESS(COMMAND ${OCTAVE_CONFIG} -p OCTLIBDIR
      OUTPUT_VARIABLE OCTLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

    # Make the values accessible from other CMakeLists.txt files
    # Also, this allows packagers to override the default values
    SET(FEELPP_OCT_DIR ${OCT_SITE_DIR}/feel++ CACHE PATH ".oct files from Feel++")
    SET(FEELPP_M_DIR ${M_SITE_DIR}/feel++ CACHE PATH ".m files from Feel++")

    message(STATUS "[feelpp] oct dir: ${FEELPP_OCT_DIR}" )
    message(STATUS "[feelpp] m dir: ${FEELPP_M_DIR}" )
    message(STATUS "[feelpp] include dir: ${OCTINCLUDEDIR}" )


    INCLUDE_DIRECTORIES( ${OCTINCLUDEDIR} )

    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Octave" )
  endif( OCTAVE_FOUND )
endif( FEELPP_ENABLE_OCTAVE)

#
# Gmsh
#
if(FEELPP_USE_GMSH_PACKAGE)
	FIND_PACKAGE(Gmsh)
else()
	set(GMSH_FOUND false)
endif()
if(NOT GMSH_FOUND)#Download and Instal it
  message(STATUS "[feelpp] GMSH NOT FOUND - Downloading and Installing it" )
  execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/contrib/gmsh-compile)
  message(STATUS "[feelpp] Building gmsh in ${CMAKE_BINARY_DIR}/contrib/gmsh-compile...")
  execute_process(
    COMMAND ${FEELPP_HOME_DIR}/contrib/gmsh/gmsh.sh ${CMAKE_BINARY_DIR}/contrib/gmsh ${NProcs2}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/gmsh-compile
    #      OUTPUT_QUIET
    OUTPUT_FILE "gmsh-configure"
    )

  FIND_PACKAGE(Gmsh REQUIRED)
endif()
if ( GMSH_FOUND )
  ADD_DEFINITIONS( -DFEELPP_HAS_GMSH=1 -D_FEELPP_HAS_GMSH_ -DGMSH_EXECUTABLE=${GMSH_EXECUTABLE} )
  if ( GL2PS_LIBRARY )
    if ( GL_LIBRARY AND FEELPP_ENABLE_OPENGL )
      SET(FEELPP_LIBRARIES ${GMSH_LIBRARIES} ${GL2PS_LIBRARY} ${GL_LIBRARY} ${FEELPP_LIBRARIES})
    else()
      SET(FEELPP_LIBRARIES ${GMSH_LIBRARIES} ${GL2PS_LIBRARY} ${FEELPP_LIBRARIES})
    endif()
  else()
    SET(FEELPP_LIBRARIES ${GMSH_LIBRARIES} ${FEELPP_LIBRARIES})
  endif()
  include_directories(${GMSH_INCLUDE_DIR})
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Gmsh" )
endif()
include(feelpp.module.gmsh)

#
# if Feel++ has been installed on the system
#
if ( NOT EXISTS ${CMAKE_SOURCE_DIR}/feel OR NOT EXISTS ${CMAKE_SOURCE_DIR}/contrib )
  include(feelpp.macros)
  FIND_PATH(FEELPP_INCLUDE_DIR feel/feelconfig.h  PATHS $ENV{FEELPP_DIR}/include/ PATH_SUFFIXES feel NO_DEFAULT_PATH )
  FIND_PATH(FEELPP_INCLUDE_DIR feel/feelconfig.h  PATHS /usr/include /opt/local/include PATH_SUFFIXES feel )

  #  FIND_LIBRARY(FEELPP_GFLAGS_LIBRARY feelpp_gflags PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  #  FIND_LIBRARY(FEELPP_GLOG_LIBRARY feelpp_glog PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  #  FIND_LIBRARY(FEELPP_CLN_LIBRARY feelpp_cln PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  if ( FEELPP_ENABLE_NLOPT )
    FIND_LIBRARY(FEELPP_NLOPT_LIBRARY feelpp_nlopt PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  endif()
  FIND_LIBRARY(FEELPP_GINAC_LIBRARY feelpp_ginac PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  FIND_LIBRARY(FEELPP_LIBRARY feelpp PATHS $ENV{FEELPP_DIR}/lib NO_DEFAULT_PATH)
  FIND_LIBRARY(FEELPP_LIBRARY feelpp )

  INCLUDE_DIRECTORIES ( ${FEELPP_INCLUDE_DIR} ${FEELPP_INCLUDE_DIR}/feel )
  FIND_PACKAGE_HANDLE_STANDARD_ARGS (Feel DEFAULT_MSG
    FEELPP_INCLUDE_DIR  FEELPP_LIBRARY
    )

  FIND_PATH( FEELPP_DATADIR cmake/modules/FindFeel++.cmake
    PATH $ENV{FEELPP_DIR}/share/feel /usr/share/feel /usr/local/share/feel )

  message(STATUS "[feelpp] Feel++ includes: ${FEELPP_INCLUDE_DIR}")
  message(STATUS "[feelpp] Feel++ library: ${FEELPP_LIBRARY}")
  message(STATUS "[feelpp] Feel++ data: ${FEELPP_DATADIR}")
  
  MARK_AS_ADVANCED(
    FEELPP_INCLUDE_DIR
    FEELPP_LIBRARY
    )
  SET(FEELPP_LIBRARIES ${FEELPP_LIBRARY} ${FEELPP_GINAC_LIBRARY} ${FEELPP_NLOPT_LIBRARY}  ${FEELPP_LIBRARIES})
else()
  set(FEELPP_LIBRARY feelpp) 
  INCLUDE_DIRECTORIES (
    ${FEELPP_BUILD_DIR}/
    ${FEELPP_SOURCE_DIR}/
    ${FEELPP_SOURCE_DIR}/contrib/gmm/include
    )
endif()



LINK_DIRECTORIES(
  ${VTK_LIBRARY_DIRS}
  ${BOOST_LIBRARY_PATH}
  ${MPI_LIBRARY_PATH}
  )


MARK_AS_ADVANCED(FEELPP_LIBRARIES)
