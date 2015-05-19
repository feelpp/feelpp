cmake_minimum_required(VERSION 2.8)
project(ParMETIS)

set(GKLIB_PATH ${FEELPP_SOURCE_DIR}/contrib/metis/GKlib CACHE PATH "path to GKlib")
set(METIS_PATH ${FEELPP_SOURCE_DIR}/contrib/metis/ CACHE PATH "path to METIS")

# Symlink ./metis to wherever metis is. This allows files to be
# included from metis/libmetis/.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${METIS_PATH} metis)

# Search for MPI.
# GK commented this out as it seems to be creating problems
# include(FindMPI)
# if(NOT MPI_FOUND)
#   message(FATAL_ERROR "mpi is not found")
# endif()
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")

# Prepare libraries.
if(SHARED)
  set(ParMETIS_LIBRARY_TYPE SHARED)
else()
  set(ParMETIS_LIBRARY_TYPE STATIC)
endif()

include(${GKLIB_PATH}/GKlibSystem.cmake)

# List of paths that the compiler will search for header files.
# i.e., the -I equivalent
include_directories(include)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${GKLIB_PATH})
include_directories(${METIS_PATH}/include)

# List of directories that cmake will look for CMakeLists.txt
add_subdirectory(${METIS_PATH}/libmetis ${CMAKE_BINARY_DIR}/libmetis)
add_subdirectory(include)
add_subdirectory(libparmetis)


# This is for testing during development and is not being distributed
#add_subdirectory(test)
