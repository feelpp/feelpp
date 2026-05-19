# - Find the FFTW library
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   FFTW_FOUND               ... true if fftw is found on the system
#   FFTW_LIBRARIES           ... full path to fftw library
#   FFTW_INCLUDES            ... fftw include directory
#
# The following variables will be checked by the function
#   FFTW_USE_STATIC_LIBS    ... if true, only static libraries are found
#   FFTW_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#   FFTW_LIBRARY            ... fftw library to use
#   FFTW_INCLUDE_DIR        ... fftw include directory
#

#If environment variable FFTWDIR is specified, it has same effect as FFTW_ROOT
if( NOT FFTW_ROOT AND DEFINED ENV{FFTWDIR} )
  set( FFTW_ROOT $ENV{FFTWDIR} )
endif()

if(DEFINED ENV{SPACK_ENV} AND NOT "$ENV{SPACK_ENV}" STREQUAL "")
  list(APPEND _FFTW_PREFIX_HINTS ${CMAKE_PREFIX_PATH})
  if(DEFINED ENV{CMAKE_PREFIX_PATH} AND NOT "$ENV{CMAKE_PREFIX_PATH}" STREQUAL "")
    set(_FFTW_ENV_PREFIX_HINTS "$ENV{CMAKE_PREFIX_PATH}")
    if(NOT WIN32)
      string(REPLACE ":" ";" _FFTW_ENV_PREFIX_HINTS "${_FFTW_ENV_PREFIX_HINTS}")
    endif()
    list(APPEND _FFTW_PREFIX_HINTS ${_FFTW_ENV_PREFIX_HINTS})
  endif()
  list(FILTER _FFTW_PREFIX_HINTS EXCLUDE REGEX "^$")
  list(REMOVE_DUPLICATES _FFTW_PREFIX_HINTS)

  set(_FFTW_REAL_PREFIX_HINTS)
  foreach(_FFTW_PREFIX_HINT IN LISTS _FFTW_PREFIX_HINTS)
    if(EXISTS "${_FFTW_PREFIX_HINT}")
      get_filename_component(_FFTW_REAL_PREFIX_HINT "${_FFTW_PREFIX_HINT}" REALPATH)
      list(APPEND _FFTW_REAL_PREFIX_HINTS "${_FFTW_REAL_PREFIX_HINT}")
    endif()
  endforeach()

  foreach(_FFTW_CACHE_VAR FFTW_LIB FFTW_MPI_LIB FFTWF_LIB FFTWL_LIB FFTW_INCLUDES)
    if(DEFINED ${_FFTW_CACHE_VAR} AND NOT "${${_FFTW_CACHE_VAR}}" MATCHES "-NOTFOUND$")
      get_filename_component(_FFTW_REAL_CACHE_VALUE "${${_FFTW_CACHE_VAR}}" REALPATH)
      set(_FFTW_CACHE_IN_ACTIVE_PREFIX FALSE)
      foreach(_FFTW_REAL_PREFIX_HINT IN LISTS _FFTW_REAL_PREFIX_HINTS)
        if(_FFTW_REAL_CACHE_VALUE MATCHES "^${_FFTW_REAL_PREFIX_HINT}(/|$)")
          set(_FFTW_CACHE_IN_ACTIVE_PREFIX TRUE)
        endif()
      endforeach()
      if(NOT _FFTW_CACHE_IN_ACTIVE_PREFIX)
        message(STATUS "[fftw] ignoring cached ${_FFTW_CACHE_VAR}='${${_FFTW_CACHE_VAR}}' outside active Spack prefixes")
        unset(${_FFTW_CACHE_VAR} CACHE)
        unset(${_FFTW_CACHE_VAR})
      endif()
    endif()
  endforeach()
endif()

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT FFTW_ROOT )
  pkg_check_modules( PKG_FFTW QUIET "fftw3" )
endif()

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

if( ${FFTW_USE_STATIC_LIBS} )
  set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
else()
  set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX} )
endif()

if( FFTW_ROOT )

  #find libs
  find_library(
    FFTW_LIB
    NAMES "fftw3"
    PATHS ${FFTW_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
  )

  find_library(
    FFTW_MPI_LIB
    NAMES "fftw3_mpi"
    PATHS ${FFTW_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
  )

  find_library(
    FFTWF_LIB
    NAMES "fftw3f"
    PATHS ${FFTW_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
  )

  find_library(
    FFTWL_LIB
    NAMES "fftw3l"
    PATHS ${FFTW_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
  )

  #find includes
  find_path(
    FFTW_INCLUDES
    NAMES "fftw3.h"
    PATHS ${FFTW_ROOT}
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
  )

else()

  find_library(
    FFTW_LIB
    NAMES "fftw3"
    PATHS ${_FFTW_PREFIX_HINTS} ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    PATH_SUFFIXES "lib" "lib64"
  )

  find_library(
    FFTW_MPI_LIB
    NAMES "fftw3_mpi"
    PATHS ${_FFTW_PREFIX_HINTS} ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    PATH_SUFFIXES "lib" "lib64"
  )

  find_library(
    FFTWF_LIB
    NAMES "fftw3f"
    PATHS ${_FFTW_PREFIX_HINTS} ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    PATH_SUFFIXES "lib" "lib64"
  )


  find_library(
    FFTWL_LIB
    NAMES "fftw3l"
    PATHS ${_FFTW_PREFIX_HINTS} ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    PATH_SUFFIXES "lib" "lib64"
  )

  find_path(
    FFTW_INCLUDES
    NAMES "fftw3.h"
    PATHS ${_FFTW_PREFIX_HINTS} ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES "include"
  )

endif( FFTW_ROOT )

if(FFTW_LIB)
  set(FFTW_LIBRARIES ${FFTW_LIB})
endif()

if(FFTW_MPI_LIB)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_MPI_LIB})
endif()

if(FFTWF_LIB)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTWF_LIB})
endif()

if(FFTWL_LIB)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTWL_LIB})
endif()

MESSAGE(STATUS "FFTW_INCLUDES = ${FFTW_INCLUDES}")
MESSAGE(STATUS "FFTW_LIBRARIES = ${FFTW_LIBRARIES}")
MESSAGE(STATUS "FFTW_LIB = ${FFTW_LIB}")
MESSAGE(STATUS "FFTW_MPI_LIB = ${FFTW_MPI_LIB}")
MESSAGE(STATUS "FFTWF_LIB = ${FFTWF_LIB}")
MESSAGE(STATUS "FFTWL_LIB = ${FFTWL_LIB}")

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
                                  FFTW_INCLUDES FFTW_LIBRARIES)

mark_as_advanced(FFTW_INCLUDES FFTW_LIBRARIES FFTW_LIB FFTWF_LIB FFTWL_LIB)
