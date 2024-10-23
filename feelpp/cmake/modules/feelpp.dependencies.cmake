# - Find Feel
# This module looks for Feel (Library for the Finite Element Method) support
# it will define the following values
#  FEELPP_INCLUDE_DIR = where feel/feelcore/feel.hpp can be found
#  FEELPP_LIBRARY    = the library to link in
#  FEELPP_LIBRARIES = list of depend libraries

# define the feel++ c++ standard level, it used to be hardcoded, this way we can
# have builds to test the different standard flavors
if (NOT DEFINED FEELPP_STD_CPP )
  set(FEELPP_STD_CPP "17") # DOC STRING "define feel++ standard c++ (default c++11), values can be : 11, 14, 17, 2a")
endif()
if (NOT DEFINED FEELPP_STDLIB_CPP AND NOT APPLE)
  set(FEELPP_STDLIB_CPP "stdc++") # DOC STRING "define feel++ standard c++ library (default libstdc++), values can be : libc++ libstdc++")
elseif( NOT DEFINED FEELPP_STDLIB_CPP )
  set(FEELPP_STDLIB_CPP "c++") # DOC STRING "define feel++ standard c++ library (default libstdc++), values can be : libc++ libstdc++")
endif()
message(STATUS "[feelpp] using c++${FEELPP_STD_CPP} standard." )
message(STATUS "[feelpp] using lib${FEELPP_STDLIB_CPP} standard c++ library." )
# Check compiler
message(STATUS "[feelpp] Compiler version : ${CMAKE_CXX_COMPILER_ID}  ${CMAKE_CXX_COMPILER_VERSION}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-cpp -Wno-deprecated-declarations" )
  # require at least gcc 4.9
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
      message(WARNING "GCC version must be at least 4.9. it is currently version ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # require at least clang 3.4
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
      string(COMPARE EQUAL "${CMAKE_CXX_COMPILER_VERSION}" "" CLANG_VERSION_EMPTY)
      if(CLANG_VERSION_EMPTY)
          message(WARNING "CMake was unable to check Clang version. It will assume that the version requirements are met.")
      else()
          message(FATAL_ERROR "Clang version must be at least 3.4! we have clang ${CMAKE_CXX_COMPILER_VERSION}")
      endif()
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

  #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${FEELPP_STD_CPP}" )
  #set(CMAKE_CXX_STANDARD ${FEELPP_STD_CPP})
  #set(CMAKE_CXX_STANDARD_REQUIRED ON)

  if ( NOT ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel") )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=1024" )
  endif()
  if ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd3373" )
  endif()
  IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=lib${FEELPP_STDLIB_CPP}" )
    if ( "${FEELPP_STDLIB_CPP}" STREQUAL "c++" )
      find_path(CXXABI_H cxxabi.h PATH_SUFFIXES libcxxabi)
      message(STATUS "[feelpp] use ${CXXABI_H}")
      if ( CXXABI_H )
        include_directories(${CXXABI_H})
      else()
        message(ERROR "[feelpp] ${CXXABI_H}")
      endif()

    ENDIF()
  endif()
ENDIF()

# IBM XL compiler
IF( ("${CMAKE_CXX_COMPILER_ID}" MATCHES "XL") )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qlanglvl=extc1x" )
endif()

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden") # -fvisibility=hidden" )

# Sometimes relinking libraries in contrib does not work becaus of bad paths
# A discussion has been opened on the bugtracker of CMake: https://cmake.org/Bug/print_bug_page.php?bug_id=13934
# To fix this: if the executable format has not been specified, which seems to happen with Clang, we force it to ELF for Linux
# Potentially also a bug to fix on OS X if it happens (The format is MACHO on OS X)
if(NOT CMAKE_EXECUTABLE_FORMAT)
    if(CMAKE_SYSTEM_NAME MATCHES "Linux")
        message(STATUS "CMAKE_EXECUTABLE_FORMAT is not set. Setting to ELF on current system (Linux).")
        set(CMAKE_EXECUTABLE_FORMAT "ELF")
    else()
        message(WARNING "CMAKE_EXECUTABLE_FORMAT is not set, you might end up with relinking errors with contrib libraries.")
    endif()
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
INCLUDE(CMakeDependentOption)

OPTION(FEELPP_ENABLE_SYSTEM_EIGEN3 "enable system eigen3 support" OFF)



OPTION(FEELPP_ENABLE_MOVE_SEMANTICS "enable move semantics(elision)" ON )
OPTION(FEELPP_ENABLE_INSTANTIATION_MODE "Instantiation mode" ON )
OPTION(FEELPP_ENABLE_MPI_MODE "Instantiation mode" ON )
OPTION(FEELPP_ENABLE_SCHED_SLURM "Enable Feel++ slurm submission scripts generation" OFF)
OPTION(FEELPP_ENABLE_SCHED_CCC "Enable Feel++ tgcc/ccc submission scripts generation" OFF)
OPTION(FEELPP_ENABLE_SCHED_LOADLEVELER "Enable Feel++ ibm(supermuc) submission scripts generation" OFF)
OPTION(FEELPP_ENABLE_TBB "enable feel++ TBB support" OFF)
OPTION(FEELPP_ENABLE_TRILINOS "enable feel++ Trilinos support" OFF)
OPTION(FEELPP_ENABLE_EXODUS "enable feel++ Exodus support" OFF)
#if ( APPLE )
  #OPTION(FEELPP_ENABLE_OPENTURNS "enable feel++ OpenTURNS support" OFF)
#else()
  #OPTION(FEELPP_ENABLE_OPENTURNS "enable feel++ OpenTURNS support" ON)
#endif()
OPTION(FEELPP_ENABLE_OCTAVE "Enable Feel++/Octave interface" OFF)


OPTION(FEELPP_ENABLE_OPENGL "enable feel++ OpenGL support" ON)
OPTION(FEELPP_DISABLE_EIGEN_ALIGNMENT "disable alignement (hence vectorization) in Eigen" OFF)

if(FEELPP_MINIMAL_BUILD)
  set( FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION OFF)
else()
  set( FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION ON)
endif()

function(get_linux_lsb_release_information)
    find_program(LSB_RELEASE_EXEC lsb_release)
    if(NOT LSB_RELEASE_EXEC)
        message(WARNING "Could not detect lsb_release executable, can not gather required information")
    else()
      execute_process(COMMAND "${LSB_RELEASE_EXEC}" --short --id OUTPUT_VARIABLE LSB_RELEASE_ID_SHORT OUTPUT_STRIP_TRAILING_WHITESPACE)
      execute_process(COMMAND "${LSB_RELEASE_EXEC}" --short --release OUTPUT_VARIABLE LSB_RELEASE_VERSION_SHORT OUTPUT_STRIP_TRAILING_WHITESPACE)
      execute_process(COMMAND "${LSB_RELEASE_EXEC}" --short --codename OUTPUT_VARIABLE LSB_RELEASE_CODENAME_SHORT OUTPUT_STRIP_TRAILING_WHITESPACE)

      set(LSB_RELEASE_ID_SHORT "${LSB_RELEASE_ID_SHORT}" PARENT_SCOPE)
      set(LSB_RELEASE_VERSION_SHORT "${LSB_RELEASE_VERSION_SHORT}" PARENT_SCOPE)
      set(LSB_RELEASE_CODENAME_SHORT "${LSB_RELEASE_CODENAME_SHORT}" PARENT_SCOPE)
    endif()
endfunction()

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  get_linux_lsb_release_information()
  message(STATUS "Linux ${LSB_RELEASE_ID_SHORT} ${LSB_RELEASE_VERSION_SHORT} ${LSB_RELEASE_CODENAME_SHORT}")
endif()

find_package(PkgConfig REQUIRED)

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
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET( FEELPP_INSTANTIATION_MODE 1 )
  else()
    SET( FEELPP_INSTANTIATION_MODE 1 )
  endif()
ENDIF()
SET(FEELPP_MESH_MAX_ORDER "2" CACHE STRING "maximum geometrical order in templates to instantiate up to 5 in 2D and 4 in 3D" )

# enable host specific
include(feelpp.host)
if ( ${CMAKE_VERSION} VERSION_GREATER 3.1 )
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  if ( THREADS_HAVE_PTHREAD_ARG )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread" )
  endif()
  if ( CMAKE_THREAD_LIBS_INIT )
    list(APPEND FEELPP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
  endif()
else()
  find_package(Threads REQUIRED)
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread" )
  SET(FEELPP_LIBRARIES ${FEELPP_LIBRARIES} pthread)
endif()

if ( FEELPP_ENABLE_TBB )
  FIND_PACKAGE(TBB)
  IF ( TBB_FOUND )
    INCLUDE_DIRECTORIES( ${TBB_INCLUDE_DIR} )
    SET(FEELPP_LIBRARIES ${TBB_LIBRARIES} ${FEELPP_LIBRARIES})
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

# Disable searching for MPI-2 C++ bindings
set(MPI_CXX_SKIP_MPICXX TRUE)
FIND_PACKAGE(MPI REQUIRED)
IF ( MPI_FOUND )
  #SET(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")
  SET( FEELPP_HAS_MPI 1 )
  SET( FEELPP_HAS_MPI_H 1 )
  #ADD_DEFINITIONS( -DFEELPP_HAS_MPI=1 -DFEELPP_HAS_MPI_H=1 )
  SET(FEELPP_BOOST_MPI mpi)
  SET(FEELPP_LIBRARIES ${MPI_LIBRARIES} ${FEELPP_LIBRARIES})
  #INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Mpi" )

  # Check for MPI IO Support

  #TRY_COMPILE(MPIIO_SUCCESS ${CMAKE_CURRENT_BINARY_DIR}/tryCompileMPIIO
  #${CMAKE_CURRENT_SOURCE_DIR}/cmake/codes/try-mpiio.cpp
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
  # This minimal sample test for the 2.2 Types (produces garbage if tested)
  CHECK_CXX_SOURCE_COMPILES(
      "
      #include <stdint.h>
      #include <mpi.h>
      #define SIZE 64
      int main(int argc, char** argv)
      {
          int32_t buf32[SIZE];
          int64_t buf64[SIZE];
          MPI_File file;
          MPI_Status status;

          MPI_Init(&argc, &argv);
          MPI_File_open( MPI_COMM_WORLD, \"output.bin\", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file );
          MPI_File_write( file, buf32, SIZE, MPI_INT32_T, &status );
          MPI_File_write( file, buf64, SIZE, MPI_INT64_T, &status );
          MPI_File_close( &file );
          MPI_Finalize();
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

  # Find mpi.h file in an attempt to detect flavour
  find_path(__MPI_H_PATH mpi.h
      PATH ${MPI_INCLUDE_PATH}
  )
  # If we find the header we attempt to grep lines in the files
  # To detect the flavour
  if(__MPI_H_PATH)
      # Attempt to grep OPEN_MPI in the file
      file(STRINGS ${__MPI_H_PATH}/mpi.h relines REGEX "OPEN_MPI")
      if(relines)
          # Put OpenMPI specific code here

          # Automatically add oarsh to the options when using OpenMPI
          # See: https://oar.readthedocs.org/en/2.5/user/usecases.html#using-mpi-with-oarsh
          if(FEELPP_ENABLE_SCHED_OAR)
              set(MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} -mca plm_rsh_agent \"oarsh\")
              message(STATUS "[feelpp] OAR detected - Automatically adding transport option to MPIEXEC_PREFLAGS (-mca plm_rsh_agent \"oarsh\")")
          endif()
      endif()
  endif()
  unset(__MPI_H_PATH)

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

if ( 0 )
find_package(GMP)
if ( GMP_FOUND )
  SET(FEELPP_LIBRARIES  ${GMP_LIBRARIES} ${FEELPP_LIBRARIES})
  message(STATUS "[feelpp] GMP: ${GMP_LIBRARIES}" )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Gmp" )
endif()

find_package(GMM)
if ( GMM_FOUND )
  message(STATUS "[feelpp] GMM includes: ${GMM_INCLUDE_DIR}" )
  include_directories(${GMM_INCLUDE_DIR})
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Gmm" )
endif()
endif()
#
# Intel MKL
#
if (NOT DEFINED FEELPP_ENABLE_MKL ) # not necessary with cmake >= 3.13
  OPTION( FEELPP_ENABLE_MKL "Enable the Intel MKL library" OFF )
endif()
if ( FEELPP_ENABLE_MKL )
  find_package(MKL)
  if ( MKL_FOUND )

    message(STATUS "[feelpp] MKL Includes: ${MKL_INCLUDE_DIRS}")
    message(STATUS "[feelpp] MKL Libraries: ${MKL_LIBRARIES}")
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Intel(MKL)" )
 #  else( MKL_FOUND )
#     #
#     # Blas and Lapack
#     #
#     if (APPLE)
#       # FIND_LIBRARY(ATLAS_LIBRARY
#       #   NAMES
#       #   atlas
#       #   PATHS
#       #   /opt/local/lib/lib
#       #   NO_DEFAULT_PATH
#       #   )
#       # message(STATUS "[feelpp] ATLAS: ${ATLAS_LIBRARY}" )
#       # IF( ATLAS_LIBRARY )
#       #   SET(FEELPP_LIBRARIES ${ATLAS_LIBRARY} ${FEELPP_LIBRARIES})
#       # ENDIF()
#       FIND_PACKAGE(LAPACK)
#     else (APPLE)
#       FIND_PACKAGE(LAPACK)
#     endif (APPLE)
#     SET(FEELPP_LIBRARIES  ${LAPACK_LIBRARIES} ${FEELPP_LIBRARIES})
   endif(MKL_FOUND)
endif(FEELPP_ENABLE_MKL)

# HDF5
# On debian,
# - do not install hdf5-helpers, otherwise it will pick the serial version by default
# - install only the libhdf5-openmpi-dev package
set(CHECK_H5_PARALLEL_CODE "
#include <hdf5.h>
int main() {
#ifdef H5_HAVE_PARALLEL
return 0; // Parallel support enabled
#else
return 1; // Parallel support not enabled
#endif
}")
function(CheckHDF5Parallel _result_var)
  # Write the check code to a file
  file(WRITE "${CMAKE_BINARY_DIR}/check_h5_parallel.c" "${CHECK_H5_PARALLEL_CODE}")
  #message(STATUS "-DINCLUDE_DIRECTORIES=${HDF5_INCLUDE_DIRS};${MPI_INCLUDE_PATH}")
  #message(STATUS "LINK_LIBRARIES ${HDF5_LINK_LIBRARIES}")
  # Attempt to compile the test code
  try_compile(HDF5_PARALLEL
    "${CMAKE_BINARY_DIR}"
    "${CMAKE_BINARY_DIR}/check_h5_parallel.c"
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${HDF5_INCLUDE_DIRS};${MPI_INCLUDE_PATH}"
    LINK_LIBRARIES ${HDF5_LINK_LIBRARIES}
    OUTPUT_VARIABLE COMPILE_OUTPUT
    )
  # message(STATUS "COMPILE_OUTPUT=${COMPILE_OUTPUT}")
  # Set the result variable to TRUE if the code compiled successfully (indicating parallel support)
  if(HDF5_PARALLEL)
    set(${_result_var} TRUE PARENT_SCOPE)
  else()
    set(${_result_var} FALSE PARENT_SCOPE)
  endif()
endfunction()

option( FEELPP_ENABLE_HDF5 "Enable HDF5 Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if ( FEELPP_ENABLE_HDF5 )

  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 COMPONENTS C)
  if( NOT HDF5_FOUND )
    pkg_check_modules(HDF5 hdf5 IMPORTED_TARGET)
    if ( HDF5_FOUND )
      # message(STATUS "[feelpp] (pkgconfig) HDF5 - Headers ${HDF5_INCLUDE_DIRS}" )
      # message(STATUS "[feelpp] (pkgconfig) HDF5 - Libraries ${HDF5_LIBRARIES}" )
      # message(STATUS "[feelpp] (pkgconfig) HDF5 - Link Libraries ${HDF5_LINK_LIBRARIES}")
      # set(FEELPP_HAS_HDF5 1)
      # set(FEELPP_LIBRARIES ${HDF5_LIBRARIES} ${FEELPP_LIBRARIES})
      set(HDF5_PKGCONFIG TRUE)
    endif()
  endif()

  if( HDF5_FOUND )

    if( NOT HDF5_IS_PARALLEL )
      CheckHDF5Parallel(HDF5_IS_PARALLEL)
    endif()
    if(HDF5_IS_PARALLEL)
      message(STATUS "[feelpp] HDF5 compiled with parallel support.")
    else()
      MESSAGE(FATAL_ERROR "[feelpp] HDF5 has been found but is not parallel, HDF5 is not enabled in Feel++")
    endif()

    message(STATUS "[feelpp] HDF5 - Headers ${HDF5_INCLUDE_DIRS}" )
    message(STATUS "[feelpp] HDF5 - Libraries ${HDF5_LIBRARIES}" )
    set(FEELPP_LIBRARIES ${HDF5_LIBRARIES} ${FEELPP_LIBRARIES})
    set(FEELPP_HAS_HDF5 1)
    # check HDF5 version
    if(NOT HDF5_VERSION)
      set( HDF5_VERSION "" )
      foreach( _dir IN LISTS HDF5_INCLUDE_DIRS )
        foreach(_hdr "${_dir}/H5pubconf.h" "${_dir}/H5pubconf-64.h" "${_dir}/H5pubconf-32.h")
          if( EXISTS "${_hdr}" )
            #MESSAGE(STATUS "_hdr=${_hdr}")
            file( STRINGS "${_hdr}"
              HDF5_VERSION_DEFINE
              REGEX "^[ \t]*#[ \t]*define[ \t]+H5_VERSION[ \t]+" )
            #MESSAGE(STATUS "HDF5_VERSION_DEFINE=${HDF5_VERSION_DEFINE}")
            if( "${HDF5_VERSION_DEFINE}" MATCHES
                "H5_VERSION[ \t]+\"([0-9]+\\.[0-9]+\\.[0-9]+)(-patch([0-9]+))?\"" )
              set( HDF5_VERSION "${CMAKE_MATCH_1}" )
              if( CMAKE_MATCH_3 )
                set( HDF5_VERSION ${HDF5_VERSION}.${CMAKE_MATCH_3})
              endif()
            endif()
            #MESSAGE(STATUS "HDF5_VERSION=${HDF5_VERSION}")
            unset(HDF5_VERSION_DEFINE)
          endif()
        endforeach()
      endforeach()
    endif(NOT HDF5_VERSION)

    STRING (REGEX MATCHALL "[0-9]+" _versionComponents "${HDF5_VERSION}")
    #MESSAGE(STATUS "HDF5_VERSION=${HDF5_VERSION}")
    LIST(GET _versionComponents 0 HDF_VERSION_MAJOR_REF)
    LIST(GET _versionComponents 1 HDF_VERSION_MINOR_REF)
    LIST(GET _versionComponents 2 HDF_VERSION_RELEASE_REF)
    SET(HDF_VERSION_REF "${HDF5_VERSION}")

    # IF (NOT HDF_VERSION_MAJOR_REF EQUAL 1 OR NOT HDF_VERSION_MINOR_REF EQUAL 8)
    MESSAGE(STATUS "[feelpp] HDF5 version is ${HDF_VERSION_REF}")
    # ENDIF()
    set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} HDF5/${HDF_VERSION_REF}" )

  else(HDF5_FOUND)
    MESSAGE(FATAL_ERROR "[feelpp] no HDF5 found")
  endif( HDF5_FOUND)
endif(FEELPP_ENABLE_HDF5)


# XDMF
option( FEELPP_ENABLE_XDMF "Enable XDMF Support" OFF )
if ( FEELPP_ENABLE_XDMF )
  find_package(XDMF)
  if (XDMF_FOUND)
    message(STATUS "[feelpp] Xdmf: includedirs : ${XDMF_INCLUDE_DIR} ")
    if ( NOT TARGET XDMF::XDMF )
      message(STATUS "build xdmf::xdmf" )
      add_library( XDMF::XDMF ALIAS Xdmf )
      set_target_properties(XDMF::XDMF PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${XDMF_INCLUDE_DIRS}"
        )
    endif()
    #INCLUDE_DIRECTORIES( ${XDMF_INCLUDE_DIRS} )
    set(FEELPP_LIBRARIES ${XDMF_LIBRARIES} ${FEELPP_LIBRARIES})
    set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} XDMF" )
    message(STATUS "Found Xdmf." )
  else()
    message(STATUS "Could not find Xdmf." )
  endif (XDMF_FOUND)
endif()

# Python libs
option( FEELPP_ENABLE_PYTHON "Enable Python Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if(FEELPP_ENABLE_PYTHON)
  #
  # Python
  #
  FIND_PACKAGE(Python3 COMPONENTS Interpreter Development)
  if(Python3_FOUND)
    execute_process(COMMAND
      ${Python3_EXECUTABLE}
      -c "import sys; print(sys.version[0:3])"
      OUTPUT_VARIABLE PYTHON_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    message(STATUS "[feelpp] Found python version ${Python3_VERSION} (${PYTHON_VERSION})")
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} PythonInterp/${Python3_VERSION}" )
  endif()

  if ( Python3_FOUND )
    message(STATUS "[feelpp] PythonLibs: ${Python3_INCLUDE_DIRS} ${Python3_LIBRARIES}")

    # INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
    SET(FEELPP_LIBRARIES ${Python3_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} PythonLibs/${Python3_VERSION}" )
    set( FEELPP_HAS_PYTHON 1 )

    # Check that sympy is available
    include(FindPythonModules)
    find_python_module(sympy 1.1 FEELPP_SYMPY_FOUND)
    if ( FEELPP_SYMPY_FOUND )
      set( FEELPP_HAS_SYMPY 1 )
      message(STATUS "[feelpp] sympy (at least 1.1) has been found")
    else()
      message(STATUS "[feelpp] sympy (at least 1.1) has not been  found")
    endif()


    Find_Package(MPI4PY)
    if ( MPI4PY_FOUND )
      set( FEELPP_HAS_MPI4PY 1 )
      message(STATUS "[feelpp] mpi4py installed; headers: ${MPI4PY_INCLUDE_DIR}")
    else()
      message(STATUS "[feelpp] python wrapper will not be enabled")
    endif()

    Find_Package(PETSC4PY)
    if ( PETSC4PY_FOUND )
      set( FEELPP_HAS_PETSC4PY 1 )
      message(STATUS "[feelpp] petsc4py installed; headers: ${PETSC4PY_INCLUDE_DIR}")
    else()
      message(STATUS "[feelpp] petsc4py python wrapper will not be enabled")
    endif()

    find_program(FEELPP_MO2FMU mo2fmu HINTS "$ENV{HOME}/.local/bin" PATHS "$ENV{HOME}/.local/bin")
    if ( NOT FEELPP_MO2FMU-NOTFOUND )
      set(FEELPP_HAS_MO2FMU 1)
      message(STATUS "[feelpp] mo2fmu found: ${FEELPP_MO2FMU}")
    else()
      message(STATUS "[feelpp] mo2fmu not found")
    endif()
  endif()
  if (DEFINED PYTHON_SITE_PACKAGES)
    set (FEELPP_PYTHON_MODULE_PATH ${PYTHON_SITE_PACKAGES})
  else ()
    execute_process (COMMAND ${Python3_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
                      OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
                      RESULT_VARIABLE _PYTHON_pythonlib_result
                      OUTPUT_STRIP_TRAILING_WHITESPACE)

    if (_PYTHON_pythonlib_result)
      message (SEND_ERROR "Could not run ${Python3_EXECUTABLE}")
    endif ()

    get_filename_component (_ABS_PYTHON_MODULE_PATH ${_ABS_PYTHON_MODULE_PATH} ABSOLUTE)
    file (RELATIVE_PATH FEELPP_PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX} ${_ABS_PYTHON_MODULE_PATH})
  endif ()
  set (FEELPP_PYTHON${Python3_VERSION_MAJOR}_MODULE_PATH ${FEELPP_PYTHON_MODULE_PATH})
  message(STATUS "[feelpp] python module path: ${FEELPP_PYTHON_MODULE_PATH}")
endif(FEELPP_ENABLE_PYTHON)

option(FEELPP_ENABLE_PYTHON_WRAPPING "Enable Python wrapping implementation" ON)

# Boost
SET(BOOST_MIN_VERSION "1.65.0")

# Making consecutive calls to find_package for Boost to find optional components (boost_python for now)
# Making only one call to find_package and having one of the component not installed will mark Boost as not found

# First we try to find boost with the python components
if(FEELPP_ENABLE_PYTHON_WRAPPING)
    # FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} COMPONENTS python )
    # if(Boost_Python3_FOUND)
    #     set(FEELPP_HAS_BOOST_PYTHON 1)
    #     set(FEELPP_LIBRARIES ${Boost_PYTHON_LIBRARY} ${FEELPP_LIBRARIES})
    #     set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Boost-Python-Wrapping" )

    # else()
    #     message(FATAL_ERROR "[feelpp] Boost.Python was not found on your system (Required for Python Wrapping)." )
    #   endif()
endif()

# Then we try to find rest of the Boost components
if ( NOT Boost_ARCHITECTURE )
  set(Boost_ARCHITECTURE "-x64")
endif()
set(Boost_ADDITIONAL_VERSIONS "1.61" "1.62" "1.63" "1.64" "1.65" "1.66" "1.67" "1.68" "1.69" "1.70" "1.71")
set(BOOST_COMPONENTS_REQUIRED date_time filesystem system program_options unit_test_framework ${FEELPP_BOOST_MPI} regex serialization iostreams )
FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS ${BOOST_COMPONENTS_REQUIRED})
if(Boost_FOUND)
  IF(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "51")
    #add_definitions(-DBOOST_RESULT_OF_USE_TR1)
    #message(STATUS "[feelpp] added -DBOOST_RESULT_OF_USE_TR1" )
  ENDIF()
  IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    # ensures that boost.signals2 compiles with clang++ >= 3.1
    IF(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "52")
      #add_definitions(-DBOOST_NO_CXX11_VARIADIC_TEMPLATES)
      #message(STATUS "[feelpp] added -DBOOST_NO_CXX11_VARIADIC_TEMPLATES" )
    ELSE()
      #add_definitions(-DBOOST_NO_VARIADIC_TEMPLATES)
      #message(STATUS "[feelpp] added -DBOOST_NO_VARIADIC_TEMPLATES" )
    ENDIF()

    # This is meant to solves a compilation issue arising with Boost 1.56 and clang compilers (confirmed for 3.4 and 3.5)
    IF(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "55")
      #add_definitions(-DBOOST_PP_VARIADICS=0)
      #message(STATUS "[feelpp] added -DBOOST_PP_VARIADICS=0" )
    ENDIF()
  ENDIF()
  IF(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "59")
    #add_definitions(-DBOOST_OPTIONAL_USE_OLD_DEFINITION_OF_NONE=1)
    #message(STATUS "[feelpp] added -DBOOST_OPTIONAL_USE_OLD_DEFINITION_OF_NONE=1" )
  endif()
else()
  message(FATAL_ERROR "[feelpp] Please check your boost version - Should be at least ${BOOST_MIN_VERSION}")
endif()

IF ( FEELPP_ENABLE_MOVE_SEMANTICS AND Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION LESS "57" )
  SET( BOOST_UBLAS_MOVE_SEMANTICS 1 CACHE STRING "Enable Boost Ublas move semantics" FORCE )
#  ADD_DEFINITIONS( -DBOOST_UBLAS_MOVE_SEMANTICS )
ENDIF()


OPTION(BOOST_ENABLE_TEST_DYN_LINK "enable boost test with dynamic lib" ON)
MARK_AS_ADVANCED(BOOST_ENABLE_TEST_DYN_LINK)

set( BOOST_PARAMETER_MAX_ARITY 25 )
#set( BOOST_FILESYSTEM_VERSION 2)
set( BOOST_FILESYSTEM_VERSION 3)
if (BOOST_ENABLE_TEST_DYN_LINK)
#  add_definitions( -DBOOST_PARAMETER_MAX_ARITY=${BOOST_PARAMETER_MAX_ARITY} -DBOOST_TEST_DYN_LINK -DBOOST_FILESYSTEM_VERSION=${BOOST_FILESYSTEM_VERSION})
else (BOOST_ENABLE_TEST_DYN_LINK)
#  add_definitions( -DBOOST_PARAMETER_MAX_ARITY=${BOOST_PARAMETER_MAX_ARITY} -DBOOST_FILESYSTEM_VERSION=${BOOST_FILESYSTEM_VERSION})
endif (BOOST_ENABLE_TEST_DYN_LINK)

# undefined BOOST_UBLAS_TYPE_CHECK
#add_definitions(-UBOOST_UBLAS_TYPE_CHECK )
#add_definitions(-DBOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR)

# this fix an issue with boost filesystem: boost is usually no compiled with
# std=c++0x and we compile with it, this causes problems with the macro
# BOOST_SCOPED_ENUM macros whose behavior differs in both case and would
# generate different c++ codes and undefined references at link time.
# in a short future, this should not be necessary anymore
# IF(NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR NOT APPLE)
#   ADD_DEFINITIONS(-DBOOST_NO_SCOPED_ENUMS)
#   IF(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "51")
#     ADD_DEFINITIONS(-DBOOST_NO_CXX11_SCOPED_ENUMS)
#   endif()
# endif()

#INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR} ${BOOST_INCLUDE_PATH})

SET(FEELPP_LIBRARIES ${Boost_LIBRARIES} ${FEELPP_LIBRARIES})



#IF( GINAC_FOUND )
#  set( FEELPP_HAS_GINAC 1 )
#  INCLUDE_DIRECTORIES( GINAC_INCLUDE_DIRS )
#  SET(FEELPP_LIBRARIES ${GINAC_LIBRARIES} ${FEELPP_LIBRARIES})
#  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GINAC" )
#ENDIF()

#add_definitions(-DHAVE_LIBDL)

option( FEELPP_ENABLE_FFTW "Enable fftw Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if(FEELPP_ENABLE_FFTW)
  find_package(FFTW)
  if( FFTW_FOUND )
    set(FEELPP_HAS_FFTW 1)
    #INCLUDE_DIRECTORIES( ${FFTW_INCLUDES} )
    set(FEELPP_LIBRARIES ${FFTW_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} fftw" )
  endif()
endif()




#
# HARTS
#
OPTION( FEELPP_ENABLE_HARTS "Enable Harts (Runtime parallelization system)" OFF )
if ( FEELPP_ENABLE_HARTS )
  FIND_PACKAGE( HARTS )
  if( HARTS_FOUND )
    #SET(CMAKE_REQUIRED_INCLUDES ${HARTS_INCLUDES} ${CMAKE_REQUIRED_INCLUDES})
    #INCLUDE_DIRECTORIES( ${HARTS_INCLUDES} )
    ADD_DEFINITIONS(${HARTS_DEFINITIONS})
    SET(FEELPP_LIBRARIES ${HARTS_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(FEELPP_HAS_HARTS 1)

    OPTION( FEELPP_ENABLE_HARTS_DEBUG "Enable Harts Debugging" OFF )
    if ( FEELPP_ENABLE_HARTS_DEBUG )
        SET(FEELPP_HARTS_DEBUG 1)
        SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} HARTS(Debug)" )
    else()
        SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} HARTS" )
    endif()
  endif()
endif()


if ( FEELPP_ENABLE_EXODUS )
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/contrib/exodus-5.24/exodus/cbind/include/)
  add_subdirectory(contrib/exodus-5.24/exodus)
  #add_subdirectory(contrib/exodus-5.24/nemesis)
  set(FEELPP_HAS_EXODUS 1)
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Exodus" )
endif()

#
# Ann
#

option( FEELPP_ENABLE_ANN "Enable ANN Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if(FEELPP_ENABLE_ANN)
FIND_PACKAGE(ANN)
  if ( ANN_FOUND )
    set(FEELPP_HAS_ANN 1)
    set(FEELPP_HAS_ANN_H 1)
    #INCLUDE_DIRECTORIES( ${ANN_INCLUDE_DIR} )
    SET(FEELPP_LIBRARIES ${ANN_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} ANN" )
  endif()
endif()

#
# GLPK
#
option( FEELPP_ENABLE_GLPK "Enable GLPK Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if ( FEELPP_ENABLE_GLPK )
  FIND_PACKAGE(GLPK)
  if ( GLPK_FOUND )
    set(FEELPP_HAS_GLPK 1)
    set(FEELPP_HAS_GLPK_H 1)
    #INCLUDE_DIRECTORIES( ${GLPK_INCLUDE_DIR} )
    SET(FEELPP_LIBRARIES ${GLPK_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GLPK" )
  endif()
endif()

# google perf tools
option(FEELPP_ENABLE_GOOGLEPERFTOOLS "Enable Google Perf Tools (tcmalloc, stracktrace and profiler)" OFF)
if ( FEELPP_ENABLE_GOOGLEPERFTOOLS )
  find_package(GooglePerfTools)
  if ( GOOGLE_PERFTOOLS_FOUND )
    set(FEELPP_HAS_GPERFTOOLS 1 )
    message(STATUS "[feelpp] Google PerfTools: ${TCMALLOC_LIBRARIES} ${STACKTRACE_LIBRARIES} ${PROFILER_LIBRARIES}")
    #INCLUDE_DIRECTORIES(${GOOGLE_PERFTOOLS_INCLUDE_DIR})
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
#find_package(GFLAGS REQUIRED)

#INCLUDE_DIRECTORIES( ${GFLAGS_INCLUDE_DIR} )

# set(_paths "")
# set(_names "")
# feelpp_split_libs(${GFLAGS_LIBRARIES} _names _paths)
# SET(FEELPP_LIBRARIES ${_names} ${FEELPP_LIBRARIES})
# link_directories(${_paths})
# unset(_paths)
# unset(_names)

# if ( ${GFLAGS_INCLUDE_DIR} MATCHES "/contrib/" )
#   SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GFLAGS/Contrib" )
# else()
#   SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GFLAGS/System" )
# endif()

# # google glog
# find_package(GLOG REQUIRED)

# INCLUDE_DIRECTORIES( ${GLOG_INCLUDE_DIR} )

# set(_paths "")
# set(_names "")
# feelpp_split_libs(${GLOG_LIBRARIES} _names _paths)
# SET(FEELPP_LIBRARIES ${_names} ${FEELPP_LIBRARIES})
# link_directories(${_paths})
# unset(_paths)
# unset(_names)

# if ( ${GLOG_INCLUDE_DIR} MATCHES "/contrib/" )
#   SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GLOG/Contrib" )
# else()
#   SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GLOG/System" )
# endif()


# xml
option( FEELPP_ENABLE_LIBXML2 "Enable libxml2 Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if(FEELPP_ENABLE_LIBXML2)
  find_package(LibXml2 2.6.27)
  if ( LIBXML2_FOUND )
    message(STATUS "[feelpp] LibXml2: ${LIBXML2_INCLUDE_DIR} ${LIBXML2_LIBRARIES}")
    #SET(CMAKE_REQUIRED_INCLUDES "${LIBXML2_INCLUDE_DIR};${CMAKE_REQUIRED_INCLUDES}")
    #INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
    SET(FEELPP_LIBRARIES ${LIBXML2_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} LibXml2" )
    set( FEELPP_HAS_LIBXML2 1 )
  endif()
endif()

#
# GPU support
#
include(feelpp.module.gpu)

#
# Kokkos
#
include(feelpp.module.kokkos)
#
# Petsc
#
option( FEELPP_ENABLE_PETSC "Enable PETSc Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if(FEELPP_ENABLE_PETSC)
  include(feelpp.module.petsc)
endif()




if ( NOT FEELPP_HAS_SCOTCH )
  option( FEELPP_ENABLE_SCOTCH "Enable Scotch Support" OFF )
  if(FEELPP_ENABLE_SCOTCH)
    FIND_PACKAGE(Scotch)
    IF( SCOTCH_FOUND )
      message(STATUS "[feelpp] SCOTCH: ${SCOTCH_LIBRARIES}" )
      SET(FEELPP_LIBRARIES ${SCOTCH_LIBRARIES} ${FEELPP_LIBRARIES})
      set(FEELPP_HAS_SCOTCH 1)
    ENDIF()
  endif()
endif()

if ( NOT FEELPP_HAS_ML )
  option( FEELPP_ENABLE_ML "Enable ML Support" OFF )
  if(FEELPP_ENABLE_ML)
    find_package(ML)
    message(STATUS "[feelpp] ML: ${ML_LIBRARY}" )
    IF ( ML_FOUND )
      SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} ML" )
      #INCLUDE_DIRECTORIES(${ML_INCLUDE_DIR})
      SET(FEELPP_LIBRARIES ${ML_LIBRARY} ${FEELPP_LIBRARIES})
      set(FEELPP_HAS_ML 1)
    ENDIF()
  endif()
endif()

if ( NOT GFORTRAN_LIBRARY )
  FIND_LIBRARY(GFORTRAN_LIBRARY
    NAMES
    gfortran
    PATHS
    $ENV{LD_LIBRARY_PATH}
    /opt/local/lib
    /usr/lib/gcc/x86_64-linux-gnu/
    PATH_SUFFIXES
    gcc6 gcc5 gcc49 gcc48 gcc47 gcc46 gcc45 gcc44 4.7 4.6 4.5 4.4
    )
  message(STATUS "[feelpp] gfortran lib: ${GFORTRAN_LIBRARY} ")
  if ( GFORTRAN_LIBRARY )
    set( FEELPP_LIBRARIES ${GFORTRAN_LIBRARY} ${FEELPP_LIBRARIES})
  endif()
endif()

if ( NOT FEELPP_HAS_MUMPS )
  option( FEELPP_ENABLE_MUMPS "Enable MUMPS Support" OFF )
  if( FEELPP_ENABLE_MUMPS)
    FIND_PACKAGE(MUMPS)
    if ( GFORTRAN_LIBRARY AND MUMPS_FOUND )
      set( FEELPP_HAS_MUMPS 1 )
      set( FEELPP_LIBRARIES ${MUMPS_LIBRARIES} ${FEELPP_LIBRARIES} )
    endif()
  endif()
endif()


if (NOT FEELPP_HAS_SUITESPARSE)
  option( FEELPP_ENABLE_SUITESPARSE "Enable SuiteSparse Support" OFF )
  if(FEELPP_ENABLE_SUITESPARSE)
    FIND_LIBRARY(SUITESPARSECONFIG_LIBRARY
      NAMES
      suitesparseconfig
      PATHS
      $ENV{PETSC_DIR}/lib
      $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
      $ENV{SUITESPARSE_DIR}/lib
      )
    IF ( SUITESPARSECONFIG_LIBRARY )
      message(STATUS "[feelpp] SuiteSparseConfig: ${SUITESPARSECONFIG_LIBRARY}" )
      SET(FEELPP_LIBRARIES  ${SUITESPARSECONFIG_LIBRARY} ${FEELPP_LIBRARIES})
      set(FEELPP_HAS_SUITESPARSE 1)
    endif()
  endif()
endif()

if (NOT FEELPP_HAS_AMD_LIB)
  option( FEELPP_ENABLE_AMD "Enable AMD Library Support" OFF )
  if(FEELPP_ENABLE_AMD)
    FIND_LIBRARY(AMD_LIBRARY
      NAMES
      amd
      PATHS
      $ENV{PETSC_DIR}/lib
      $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
      $ENV{SUITESPARSE_DIR}/lib
      )
    IF ( AMD_LIBRARY )
      message(STATUS "[feelpp] Amd: ${AMD_LIBRARY}" )
      SET(FEELPP_LIBRARIES  ${AMD_LIBRARY} ${FEELPP_LIBRARIES})
      set(FEELPP_HAS_AMD_LIB 1)
    endif()
  endif()
endif()

if (NOT FEELPP_HAS_COLAMD_LIB)
  option( FEELPP_ENABLE_COLAMD "Enable COLAMD Library Support" OFF )
  if(FEELPP_ENABLE_COLAMD)
    FIND_LIBRARY(COLAMD_LIBRARY
      NAMES
      colamd
      PATHS
      $ENV{PETSC_DIR}/lib
      $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
      $ENV{SUITESPARSE_DIR}/lib
      )
    IF ( COLAMD_LIBRARY )
      message(STATUS "[feelpp] ColAmd: ${COLAMD_LIBRARY}" )
      SET(FEELPP_LIBRARIES  ${COLAMD_LIBRARY} ${FEELPP_LIBRARIES})
      set(FEELPP_HAS_COLAMD_LIB 1)
    endif()
  endif()
endif()

if (NOT FEELPP_HAS_CHOLMOD_LIB)
  option( FEELPP_ENABLE_CHOLMOD "Enable CHOLMOD Library Support" OFF )
  if(FEELPP_ENABLE_CHOLMOD)
    FIND_LIBRARY(CHOLMOD_LIBRARY
      NAMES
      cholmod
      PATHS
      $ENV{PETSC_DIR}/lib
      $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
      $ENV{SUITESPARSE_DIR}/lib
      )
    if ( CHOLMOD_LIBRARY )
      message(STATUS "[feelpp] Cholmod: ${CHOLMOD_LIBRARY}" )
      SET(FEELPP_LIBRARIES ${UMFPACK_LIBRARY} ${CHOLMOD_LIBRARY} ${FEELPP_LIBRARIES})
       set(FEELPP_HAS_CHOLMOD_LIB 1)
    endif()
  endif()
endif()

if (NOT FEELPP_HAS_UMFPACK_LIB)
  option( FEELPP_ENABLE_UMFPACK "Enable UMFPACK Library Support" OFF )
  if(FEELPP_ENABLE_UMFPACK)
    FIND_LIBRARY(UMFPACK_LIBRARY
      NAMES
      umfpack
      PATHS
      $ENV{PETSC_DIR}/lib
      $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
      $ENV{SUITESPARSE_DIR}/lib
      )
    if ( UMFPACK_LIBRARY )
      message(STATUS "[feelpp] Umfpack: ${UMFPACK_LIBRARY}" )
      SET(FEELPP_LIBRARIES ${UMFPACK_LIBRARY} ${FEELPP_LIBRARIES})
      set(FEELPP_HAS_UMFPACK_LIB 1)
    endif()
  endif()
endif()

if (NOT FEELPP_HAS_YAML)
  option( FEELPP_ENABLE_YAML "Enable YAML Library Support" OFF )
  if ( FEELPP_ENABLE_YAML )
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
       message(STATUS "[feelpp] YAML: ${YAML_LIBRARY}" )
      SET(FEELPP_LIBRARIES ${YAML_LIBRARY} ${FEELPP_LIBRARIES})
      set(FEELPP_HAS_YAML 1)
    endif()
  endif()
endif()


#
# parpack
#
option( FEELPP_ENABLE_PARPACK "Enable ParPack Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if(FEELPP_ENABLE_PARPACK)
  FIND_LIBRARY(PARPACK_LIBRARY NAMES parpack)
  if (PARPACK_LIBRARY)
    SET(PARPACK_LIBRARIES ${PARPACK_LIBRARY})
    SET(FEELPP_LIBRARIES ${PARPACK_LIBRARIES} ${FEELPP_LIBRARIES})
  endif()
  MARK_AS_ADVANCED( PARPACK_LIBRARY )
endif()


#
# SLEPc
#
option( FEELPP_ENABLE_SLEPC "Enable SLEPc Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if (FEELPP_ENABLE_SLEPC)
  FIND_PACKAGE( SLEPc )
  if ( SLEPC_FOUND )
    #SET(CMAKE_REQUIRED_INCLUDES "${SLEPC_INCLUDES};${CMAKE_REQUIRED_INCLUDES}")
    #INCLUDE_DIRECTORIES( ${SLEPC_INCLUDE_DIR} )
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
    #INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIR})
    SET(FEELPP_LIBRARIES ${TRILINOS_LIBRARIES} ${FEELPP_LIBRARIES})
    SET(BACKEND_TRILINOS trilinos)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Trilinos" )
  endif( TRILINOS_FOUND )
endif (FEELPP_ENABLE_TRILINOS)

#
# OpenTURNS
#
# option( FEELPP_ENABLE_OPENTURNS "Enable OpenTurns Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
option( FEELPP_ENABLE_OPENTURNS "Enable OpenTurns Support" OFF)
IF ( FEELPP_ENABLE_OPENTURNS )
  FIND_PACKAGE( OpenTURNS )
  if ( OPENTURNS_FOUND )
    MESSAGE(STATUS "[feelpp] OpenTURNS Libraries: ${OpenTURNS_LIBRARIES}")
    MESSAGE(STATUS "[feelpp] OpenTURNS Headers: ${OpenTURNS_INCLUDE_DIRS}")
    #INCLUDE_DIRECTORIES(${OpenTURNS_INCLUDE_DIRS})
    #SET(FEELPP_LIBRARIES ${OpenTURNS_LIBRARIES} ${FEELPP_LIBRARIES})
    # now OpenTURNS_LIBRARIES are used in crb_add_python_module
    set(FEELPP_HAS_OPENTURNS 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} OpenTURNS" )
  endif( OPENTURNS_FOUND )
endif()

#
# VTK
#
option( FEELPP_ENABLE_VTK "Enable VTK Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if ( FEELPP_ENABLE_VTK )
    # MESSAGE("Finding VTK:")
    # MESSAGE("PARAVIEW_DIR=$ENV{PARAVIEW_DIR}")
    # MESSAGE("MACHINE_PARAVIEW_DIR=${MACHINE_PARAVIEW_DIR}")
    # # if(EXISTS $ENV{PARAVIEW_DIR})
    # #  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{PARAVIEW_DIR}/Modules/")
    # # endif()
    # MESSAGE("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")

    option( FEELPP_ENABLE_VTK_FROM_PARAVIEW "Enable VTK Support" ON ) #${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
    if ( FEELPP_ENABLE_VTK_FROM_PARAVIEW )
    # First try to find ParaView
    # FIND_PACKAGE(ParaView QUIET
    #    COMPONENTS vtkParallelMPI vtkPVCatalyst vtkPVPythonCatalyst
    #    PATHS $ENV{PARAVIEW_DIR} ${MACHINE_PARAVIEW_DIR})

    FIND_PACKAGE(ParaView QUIET NO_MODULE
      PATHS $ENV{PARAVIEW_DIR} ${MACHINE_PARAVIEW_DIR} )
    endif()

    if(ParaView_FOUND)
      if ( PARAVIEW_USE_FILE )
        message(STATUS "[ParaView] Use file: ${PARAVIEW_USE_FILE}")
        INCLUDE(${PARAVIEW_USE_FILE})
      endif()
        # trying to load a minimal vtk
        IF (TARGET vtkParallelMPI)
        FIND_PACKAGE(ParaView QUIET COMPONENTS vtkParallelMPI NO_MODULE
         PATHS $ENV{PARAVIEW_DIR} ${MACHINE_PARAVIEW_DIR} )
          message(STATUS "[ParaView] Loading vtkParallelMPI module")
        ENDIF ()
        IF (TARGET vtkPVCatalyst)
        FIND_PACKAGE(ParaView COMPONENTS vtkPVCatalyst NO_MODULE
          PATHS $ENV{PARAVIEW_DIR} ${MACHINE_PARAVIEW_DIR} )
          message(STATUS "[ParaView] Loading vtkPVCatalyst module")
        ENDIF ()
        IF (TARGET vtkPVPythonCatalyst)
        FIND_PACKAGE(ParaView COMPONENTS vtkPVPythonCatalyst NO_MODULE
          PATHS $ENV{PARAVIEW_DIR} ${MACHINE_PARAVIEW_DIR} )
          message(STATUS "[ParaView] Loading vtkPVPythonCatalyst module")
        ENDIF ()

        # Enable VTK exporter and insitu in config
        IF (TARGET vtkPVCatalyst AND TARGET vtkPVPythonCatalyst)
          set(FEELPP_VTK_INSITU_ENABLED 1)
          message(STATUS "In-situ visualisation enabled (vtkPVCatalyst, vtkPVPythonCatalyst).")
        ENDIF()

        # Mark VTK and ParaView as available
        set(FEELPP_HAS_VTK 1)
        set(FEELPP_HAS_PARAVIEW 1)
        set(FEELPP_PARAVIEW_DIR ${ParaView_DIR})
        # Check for version to ensure that we are able to
        # use an external communicator
        set(VTK_HAS_PARALLEL 0)
        if( VTK_MAJOR_VERSION EQUAL 6 OR VTK_MAJOR_VERSION GREATER 6 )
          set(VTK_HAS_PARALLEL 1)
          # MESSAGE("VTK_HAS_PARALLEL=${VTK_HAS_PARALLEL}")
        endif()

        #SET(CMAKE_REQUIRED_INCLUDES "${VTK_INCLUDE_DIRS};${ParaView_INCLUDE_DIRS};${CMAKE_REQUIRED_INCLUDES}")
        #INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
        #INCLUDE_DIRECTORIES(${ParaView_INCLUDE_DIRS})
        #INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
        SET(FEELPP_LIBRARIES ${ParaView_LIBRARIES} ${FEELPP_LIBRARIES})
        SET(FEELPP_LIBRARIES ${VTK_LIBRARIES} ${FEELPP_LIBRARIES})
        # Generate FEELPP_VTK_LIBRARY and FEELPP_VTK_DIRS from linker
        #feelpp_expand_target_libraries( FEELPP_VTK ${VTK_LIBRARIES} )
        #set(FEELPP_LINK_LIBRARIES ${VTK_LIBRARIES} ${FEELPP_LINK_LIBRARIES})
        #set( FEELPP_LIBRARIES ${FEELPP_VTK_LIBRARIES} ${FEELPP_LIBRARIES} )
        SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} ParaView/VTK" )

        message(STATUS "Found ParaView ${PARAVIEW_VERSION_FULL}/VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
    endif()

    # If ParaView was not found we try to find VTK
    if(NOT ParaView_FOUND)
        message(STATUS "ParaView could not be found or is not compatible with Feel++. Looking for VTK...")
        FIND_PACKAGE(VTK QUIET)
        if( VTK_FOUND )
            include(${VTK_USE_FILE})

            set(FEELPP_HAS_VTK 1)
            set(FEELPP_VTK_DIR ${VTK_DIR})
            MESSAGE(STATUS "[feelpp] Found VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")# ${VTK_LIBRARIES}")

            # Check for MPI support in VTK
            set(VTK_HAS_PARALLEL 0)
            # Prior to VTK version 6, VTK_KITS was used
            if( VTK_MAJOR_VERSION LESS 6 )
                #message("Available VTK KITS: ${VTK_KITS}")
                list(FIND VTK_KITS "PARALLEL" __test_vtk_parallel)
                if( NOT ( ${__test_vtk_parallel} EQUAL -1 ) )
                    set(VTK_HAS_PARALLEL 1)
                    message(WARNING "External initialization of MPI Communicator is not available in VTK 5. VTK parallel export will be disabled.")
                else()
                    message(WARNING "MPI support for VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_BUILD_VERSION} is not activated. VTK parallel export will be disabled.")
                endif()
                unset(__test_vtk_parallel)
                # From version 6 modules ared used
            else()
                #message(FATAL_ERROR "${VTK_MODULES_ENABLED}")
                list(FIND VTK_MODULES_ENABLED "vtkParallelMPI" __test_vtk_parallel)
                if( NOT (${__test_vtk_parallel} EQUAL -1) )
                    set(VTK_HAS_PARALLEL 1)
                else()
                    message(WARNING "MPI support for VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_BUILD_VERSION} is not activated. VTK parallel export will be disabled.")
                endif()
                unset(__test_vtk_parallel)
            endif()

            #if ( NOT FEELPP_ENABLE_OPENGL )
                #SET(VTK_LIBRARIES "-lvtkRendering -lvtkGraphics -lvtkImaging  -lvtkFiltering -lvtkCommon -lvtksys" )
            #endif()
            #INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
            MARK_AS_ADVANCED( VTK_DIR )

            # # On debian for vtk6 (actually not working since vtk6 is built with hdf5 serial)
            # if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
            #   message (STATUS "Trying to find Linux distro")
            #   execute_process(
		    #     COMMAND lsb_release -i -s
		    #     RESULT_VARIABLE RC
		    #     OUTPUT_VARIABLE DIST_NAME
            #     OUTPUT_STRIP_TRAILING_WHITESPACE
            #     ERROR_STRIP_TRAILING_WHITESPACE
		    #     )
	        #   if (RC EQUAL 0)
		    #     execute_process(
		    #       COMMAND lsb_release -c -s
		    #       RESULT_VARIABLE RC
		    #       OUTPUT_VARIABLE DIST_CODENAME
            #       OUTPUT_STRIP_TRAILING_WHITESPACE
            #       ERROR_STRIP_TRAILING_WHITESPACE
		    #       )
		    #     message (STATUS "Distribution: ${DIST_NAME} ${DIST_CODENAME}")
	        #   endif ()

            #   set(DebianDistros  jessie stretch sid)
            #   if (DIST_NAME STREQUAL "Debian")
            #     list (FIND DebianDistros  ${DIST_CODENAME} _index)
            #     if (${_index} GREATER -1)
            #       find_package(Qt5Widgets REQUIRED)
            #       MESSAGE(STATUS "add ${Qt5Widgets_LIBRARIES} to FEELPP_LIBRARIES")
            #     endif()
            #   endif()
            # endif()

            # Generate FEELPP_VTK_LIBRARY and FEELPP_VTK_DIRS
            SET(FEELPP_LIBRARIES ${VTK_LIBRARIES} ${FEELPP_LIBRARIES})
            #feelpp_expand_target_libraries( FEELPP_VTK ${VTK_LIBRARIES} )
            #SET(FEELPP_LIBRARIES ${FEELPP_VTK_LIBRARIES} ${FEELPP_LIBRARIES})
            SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} VTK" )
        endif()
   endif()

    # If VTK was not found
    if(NOT ParaView_FOUND AND NOT VTK_FOUND)
        message(STATUS "Neither ParaView nor VTK were found. VTK exporter and In-situ processing not enabled")
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
    set(FEELPP_HAS_OCTAVE 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Octave" )
  endif( OCTAVE_FOUND )
endif( FEELPP_ENABLE_OCTAVE)

#
# Gmsh
#

option( FEELPP_ENABLE_GSL "Enable GSL Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )

if( FEELPP_ENABLE_GSL )
  FIND_PACKAGE(GSL)
  if ( GSL_FOUND )
    #ADD_DEFINITIONS( -DFEELPP_HAS_GSL=1 )
    SET(FEELPP_LIBRARIES ${GSL_LIBRARIES} ${FEELPP_LIBRARIES})
    #INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
    SET(FEELPP_HAS_GSL 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} GSL(${GSL_VERSION})" )
    message(STATUS "[feelpp] gsl ${GSL_VERSION} enabled" )
  endif()

endif(FEELPP_ENABLE_GSL)

#
# Gmsh
#

option( FEELPP_ENABLE_GMSH "Enable Gmsh Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if( FEELPP_ENABLE_GMSH )
  if(FEELPP_USE_GMSH_PACKAGE)
	FIND_PACKAGE(Gmsh)
  else()
	set(GMSH_FOUND false)
  endif()
  if(NOT GMSH_FOUND)#Download and Instal it
    message(STATUS "[feelpp] GMSH NOT FOUND - Downloading and Installing it" )
    execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/contrib/gmsh-compile)
    message(STATUS "[feelpp] Building gmsh in ${CMAKE_CURRENT_BINARY_DIR}/contrib/gmsh-compile...")
    execute_process(
      COMMAND ${FEELPP_HOME_DIR}/contrib/gmsh/gmsh.sh ${CMAKE_BINARY_DIR}/contrib/gmsh/ ${FEELPP_HOME_DIR}/contrib/gmsh/patches ${NProcs2} ${CMAKE_CXX_COMPILER}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/gmsh-compile
      #      OUTPUT_QUIET
      OUTPUT_FILE "gmsh-configure"
    )

    FIND_PACKAGE(Gmsh REQUIRED)
  endif()
  if ( GMSH_FOUND )
    #ADD_DEFINITIONS( -DFEELPP_HAS_GMSH=1 -D_FEELPP_HAS_GMSH_ -DGMSH_EXECUTABLE=${GMSH_EXECUTABLE} )
    if ( FEELPP_HAS_GMSH_ADAPT_H )
      # message(STATUS "Add -DFEELPP_HAS_GMSH_ADAPT_H option")
      #ADD_DEFINITIONS( -DFEELPP_HAS_GMSH_ADAPT_H )
    endif()
    if ( GL2PS_LIBRARY )
      if ( GL_LIBRARY AND FEELPP_ENABLE_OPENGL )
        SET(FEELPP_LIBRARIES ${GMSH_LIBRARIES} ${GL2PS_LIBRARY} ${GL_LIBRARY} ${FEELPP_LIBRARIES})
      else()
        SET(FEELPP_LIBRARIES ${GMSH_LIBRARIES} ${GL2PS_LIBRARY} ${FEELPP_LIBRARIES})
      endif()
    else()
      SET(FEELPP_LIBRARIES ${GMSH_LIBRARIES} ${FEELPP_LIBRARIES})
    endif()
    #INCLUDE_DIRECTORIES(${GMSH_INCLUDE_DIR})
    SET(FEELPP_HAS_GMSH 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Gmsh" )
  endif()
  include(feelpp.module.gmsh)
endif()

#
# conntrib configuration
#
#include( contrib.dependencies )

#
# Acusim
#
include(feelpp.module.altair)

# Asciidoctor
option( FEELPP_ENABLE_ASCIIDOCTOR "Enable AsciiDoctor Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if ( FEELPP_ENABLE_ASCIIDOCTOR )
  include( feelpp.adoc )
endif()

# Enable precompiled headers (PCH)
option( FEELPP_ENABLE_PCH "Enable precompiled headers (pch)" OFF )
option( FEELPP_ENABLE_PCH_APPLICATIONS "Enable precompiled headers (pch) for applications" OFF )

if( FEELPP_ENABLE_PCH )
    set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} PCH" )
endif()
if( FEELPP_ENABLE_PCH_APPLICATIONS )
    set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} PCH_Apps" )
endif()

# Enable Feel++ interpreter using cling.
option( FEELPP_ENABLE_CLING_INTERPRETER "Enable feel++ interpreter [ EXPERIMENTAL ]" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if( FEELPP_ENABLE_CLING_INTERPRETER )
  include(feelpp.module.cling.interpreter)
endif()

option( FEELPP_ENABLE_OMC "Enable OpenModelica in Feel++" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if( FEELPP_ENABLE_OMC )
  include( feelpp.module.omc )
endif()

#option( FEELPP_ENABLE_FMILIB "Enable FMILib in Feel++" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
#if( FEELPP_ENABLE_FMILIB )
#  include( feelpp.module.fmilib )
#endif()

option( FEELPP_ENABLE_LIBCURL "Enable libcurl in Feel++" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )
if ( FEELPP_ENABLE_LIBCURL )
  find_package( CURL )
  #message( "CURL_FOUND=${CURL_FOUND}" )
  #message( "CURL_INCLUDE_DIRS=${CURL_INCLUDE_DIRS}" )
  #message( "CURL_LIBRARIES=${CURL_LIBRARIES}" )
  #message( "CURL_VERSION_STRING=${CURL_VERSION_STRING}" )
  if( CURL_FOUND )
    #include_directories( ${CURL_INCLUDE_DIRS} )
    set( FEELPP_HAS_LIBCURL 1 )
    set( FEELPP_LIBRARIES ${CURL_LIBRARIES} ${FEELPP_LIBRARIES} )
    set( FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} libcurl/${CURL_VERSION_STRING}" )
  endif()
endif()

#
# if Feel++ has been installed on the system
#
if ( NOT EXISTS ${CMAKE_SOURCE_DIR}/feelpp/ OR NOT EXISTS ${CMAKE_SOURCE_DIR}/feelpp/contrib )
  include(feelpp.macros)
  FIND_PATH(FEELPP_INCLUDE_DIR feel/feelconfig.h  PATHS $ENV{FEELPP_DIR}/include/ PATH_SUFFIXES feel NO_DEFAULT_PATH )
  FIND_PATH(FEELPP_INCLUDE_DIR feel/feelconfig.h  PATHS /usr/include /opt/local/include PATH_SUFFIXES feel )

  #  FIND_LIBRARY(FEELPP_GFLAGS_LIBRARY feelpp_gflags PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  #  FIND_LIBRARY(FEELPP_GLOG_LIBRARY feelpp_glog PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  #  FIND_LIBRARY(FEELPP_CLN_LIBRARY feelpp_cln PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  find_package(CLN)
  if ( CLN_FOUND )
    INCLUDE_DIRECTORIES( ${CLN_INCLUDE_DIR} )
    SET(FEELPP_LIBRARIES ${CLN_LIBRARIES} ${FEELPP_LIBRARIES})
  endif( CLN_FOUND )

  if ( FEELPP_ENABLE_NLOPT )
    FIND_LIBRARY(FEELPP_NLOPT_LIBRARY feelpp_nlopt PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  endif()
  if ( FEELPP_ENABLE_IPOPT )
    FIND_LIBRARY(FEELPP_IPOPT_LIBRARY feelpp_ipopt feelpp_ipoptfort PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
  endif()
  if ( FEELPP_ENABLE_KWSYS )
    FIND_LIBRARY(FEELPP_KWSYS_LIBRARY feelpp_kwsys PATHS $ENV{FEELPP_DIR}/lib /usr/lib /usr/lib/feel/lib /opt/feel/lib /usr/ljk/lib )
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
SET(FEELPP_LIBRARIES ${FEELPP_LIBRARY} ${FEELPP_GINAC_LIBRARY} ${FEELPP_NLOPT_LIBRARY} ${FEELPP_IPOPT_LIBRARY} ${FEELPP_KWSYS_LIBRARY} ${FEELPP_LIBRARIES})
else()
  set(FEELPP_LIBRARY feelpp)
  SET(FEELPP_INCLUDE_DIR ${FEELPP_BUILD_DIR}/ ${FEELPP_BUILD_DIR}/feelpp ${FEELPP_SOURCE_DIR}/ ${FEELPP_SOURCE_DIR}/feelpp)
  #INCLUDE_DIRECTORIES(${FEELPP_INCLUDE_DIR})
endif()

#if ( ${FEELPP_HAS_PYTHON} AND ${FEELPP_HAS_EIGEN3} )
  #add_subdirectory( contrib/minieigen )
#endif()

# Cleaning variables.
set( varstoclean
     FEELPP_LIBRARIES )

# Do remove duplicated variable entries.
foreach( varname ${varstoclean})
    if( NOT "${${varname}}" STREQUAL "")
        list( REMOVE_DUPLICATES ${varname})
    endif()
endforeach()

MARK_AS_ADVANCED(FEELPP_LIBRARIES)
