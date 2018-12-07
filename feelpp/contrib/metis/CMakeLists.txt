cmake_minimum_required(VERSION 2.8)
project(METIS)

#SET( CMAKE_MODULE_PATH  ${FEELPP_HOME_DIR}/feelpp/contrib/metis/ )

set(C_FLAGS "-O3 -Wall -fPIC")

set(GKLIB_PATH ${CMAKE_SOURCE_DIR}/feelpp/contrib/metis/GKlib CACHE PATH "path to GKlib" )
message(STATUS "${GKLIB_PATH}")
#set(SHARED FALSE CACHE BOOL "build a shared library")

if(MSVC)
  set(METIS_INSTALL FALSE)
else()
  set(METIS_INSTALL TRUE)
endif()

# Configure libmetis library.
if(BUILD_SHARED_LIBS)
  set(METIS_LIBRARY_TYPE SHARED)
else()
  set(METIS_LIBRARY_TYPE STATIC)
endif()

include(${GKLIB_PATH}/GKlibSystem.cmake)
# Add include directories.
include_directories(${GKLIB_PATH})
include_directories(include)

# Recursively look for CMakeLists.txt in subdirs.
add_subdirectory(include)
add_subdirectory(libmetis)

