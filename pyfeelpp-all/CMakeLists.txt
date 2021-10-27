cmake_minimum_required(VERSION 3.3)
project(pyfeelpp-all)

if ( (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) OR (FEELPP_COMPONENT STREQUAL "pyfeelpp-all" ))
  find_package(Feel++ COMPONENTS Toolboxes MOR PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[pyfeelpp-all] found Feel++")
  endif()
  if(NOT FEELPP_TOOLBOXES_FOUND)
    message(FATAL_ERROR "Feel++ Toolboxes were not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[pyfeelpp-all] found Feel++ toolboxes")
  endif()
  if(NOT FEELPP_MOR_FOUND)
    message(FATAL_ERROR "Feel++ MOR were not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[pyfeelpp-all] found Feel++ MOR")
  endif()

  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "install prefix" FORCE)
    # SKIP RPATH if we set CMAKE_INSTALL_PREFIX
    SET(CMAKE_SKIP_RPATH ON)
  endif()
  
  message(STATUS "[install] ${CMAKE_INSTALL_PREFIX}")
endif()

enable_testing()

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  list(APPEND CMAKE_MODULE_PATH "$ENV{FEELPP_DIR}/share/cmake/pybind11")
  list(APPEND CMAKE_MODULE_PATH "/usr/local/share/cmake/pybind11")
  include(pybind11Tools)
endif()

add_subdirectory(pyfeelpp)
add_subdirectory(pyfeelpp-toolboxes)
add_subdirectory(pyfeelpp-mor)

