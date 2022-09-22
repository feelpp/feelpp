# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
CONFIGURE_FILE(feelconfig.h.in feelconfig.h  @ONLY)
CONFIGURE_FILE(feelinfo.h.in feelinfo.h  @ONLY)

# Define the two required variables before including
# the source code for watching a git repository.
set(PRE_CONFIGURE_FILE "feelcore/git.cpp.in")
if ( GIT_EXECUTABLE )
  set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_BINARY_DIR}/feelcore/git.cpp")
  include(GitWatcher)
elseif( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/feelcore/git.cpp")
  add_custom_target(check_git)
  set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/feelcore/git.cpp")
else()
  set(GIT_RETRIEVED_STATE "false")
  set(GIT_IS_DIRTY "false")
  add_custom_target(check_git)
  CONFIGURE_FILE(feelcore/git.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/feelcore/git.cpp @ONLY)
  set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_BINARY_DIR}/feelcore/git.cpp")
endif()
message(STATUS "[feelpp] post_configure_file: ${POST_CONFIGURE_FILE}")
SET(DIRS core timing math alg mesh partition poly opt discr ts filters material vf pde ls fit models python fmi hwsys event viewfactor)
foreach(dir ${DIRS})
  add_subdirectory(feel${dir})
endforeach()

# Add target for dynamic library
#add_library(feelpp-dynamic SHARED $<TARGET_OBJECTS:feelpp-core>)
#if(FEELPP_ENABLE_BUILD_STATIC)
#  add_library(feelpp STATIC options.cpp ${FEELPP_SRCS})
#else()
  add_library(feelpp SHARED options.cpp ${POST_CONFIGURE_FILE} ${FEELPP_SRCS})
  set_target_properties(feelpp PROPERTIES VERSION ${FEELPP_SHARED_VERSION} SOVERSION ${FEELPP_SHARED_SOVERSION})
  set_property(TARGET feelpp PROPERTY MACOSX_RPATH ON)
  set_property(TARGET feelpp PROPERTY CXX_STANDARD ${FEELPP_STD_CPP})
#endif()
add_library(Feelpp::feelpp ALIAS feelpp )  # to match exported target
add_dependencies(feelpp check_git)
set_target_properties(feelpp PROPERTIES OUTPUT_NAME "feelpp")
#set_target_properties(feelpp PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")
if ( ${CMAKE_VERSION} VERSION_GREATER 3.7.2 )
  target_compile_features(feelpp PUBLIC cxx_std_${FEELPP_STD_CPP})
endif()
set_property(TARGET feelpp PROPERTY FEELPP_STD_CPP ${FEELPP_STD_CPP} )

set_property(TARGET feelpp PROPERTY LABELS feel++)
if ( CMAKE_BUILD_TYPE MATCHES Coverage )
# target_compile_options(${project_name} INTERFACE --coverage -O0 -g)
  target_link_libraries( feelpp INTERFACE --coverage )
endif()


# put at the beginning in order to have the include ordered
target_link_libraries( feelpp PUBLIC feelpp_contrib )

target_include_directories(feelpp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/..>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:include/feelpp> )

IF ( FEELPP_ENABLE_JSON_DIAGNOSTICS )
  target_compile_definitions( feelpp PUBLIC  FEELPP_ENABLE_JSON_DIAGNOSTICS=1 )
ENDIF()

IF (FEELPP_HAS_LIBDL)
  target_link_libraries( feelpp PRIVATE ${DL_LIBRARY} )
endif()
IF (CURL_FOUND)
  target_compile_definitions( feelpp PUBLIC  FEELPP_HAS_LIBCURL )
  target_include_directories( feelpp PUBLIC ${CURL_INCLUDE_DIRS} )
  target_link_libraries( feelpp PUBLIC ${CURL_LIBRARIES} )
endif()
if ( MPI_FOUND )
  target_include_directories( feelpp PUBLIC ${MPI_INCLUDE_PATH} )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_MPI FEELPP_HAS_MPI_H )
  if ( FEELPP_HAS_MPIIO )
    target_compile_definitions( feelpp PUBLIC FEELPP_HAS_MPIIO )
  endif()
  target_link_libraries( feelpp PRIVATE ${MPI_LIBRARIES} )
  set_property( TARGET feelpp PROPERTY MPI_CXX_COMPILER ${MPI_CXX_COMPILER} )
  set_property( TARGET feelpp PROPERTY MPIEXEC ${MPIEXEC} )
  set_property( TARGET feelpp PROPERTY MPIEXEC_NUMPROC_FLAG ${MPIEXEC_NUMPROC_FLAG} )
  set_property( TARGET feelpp PROPERTY MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} )
  set_property( TARGET feelpp PROPERTY MPIEXEC_POSTFLAGS ${MPIEXEC_POSTFLAGS} )
endif()
if ( PYTHONLIBS_FOUND )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_PYTHON )
  set_property( TARGET feelpp PROPERTY FEELPP_PYTHON_MODULE_PATH ${FEELPP_PYTHON_MODULE_PATH})
  target_include_directories( feelpp PUBLIC ${PYTHON_INCLUDE_DIRS})
  target_link_libraries( feelpp PRIVATE ${PYTHON_LIBRARIES} )
endif()
if ( HDF5_FOUND AND HDF5_IS_PARALLEL )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_HDF5 )
  target_include_directories( feelpp PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries( feelpp PUBLIC ${HDF5_LIBRARIES} )
endif()
if ( XDMF_FOUND )
  target_include_directories( feelpp PRIVATE $<BUILD_INTERFACE:${XDMF_INCLUDE_DIRS}>)
  #target_link_libraries( feelpp PRIVATE ${XDMF_LIBRARIES} )
  target_link_libraries( feelpp PRIVATE XDMF::XDMF )
endif()
if ( Boost_FOUND )
  target_compile_definitions(
    feelpp
    PUBLIC
    FEELPP_HAS_BOOST
    # result_of
    BOOST_RESULT_OF_USE_TR1=1
    # Boost::none
    BOOST_OPTIONAL_USE_OLD_DEFINITION_OF_NONE=1
    # Boost::ublas
    BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR=1
    # Boost::parameter
    BOOST_PARAMETER_MAX_ARITY=${BOOST_PARAMETER_MAX_ARITY}
    # Boost::filesystem
    BOOST_FILESYSTEM_VERSION=${BOOST_FILESYSTEM_VERSION}
    )

  # Required on MacOS to circumvent error "Boost.Stacktrace requires `_Unwind_Backtrace` function"
  if( APPLE )
    target_compile_definitions(feelpp PUBLIC _GNU_SOURCE)
  endif()

  if (BOOST_ENABLE_TEST_DYN_LINK)
    target_compile_definitions( feelpp PUBLIC BOOST_TEST_DYN_LINK )
  endif (BOOST_ENABLE_TEST_DYN_LINK)

  # re-enable some std functions which have been removed in c++17 (for libcpp)
  # * std::random_shuffle : fixed in boost >= 1.69
  # * std::auto_ptr : fixed in boost >= 1.67
  if(${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION} VERSION_LESS 1.69)
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang")
      if ( "${FEELPP_STDLIB_CPP}" STREQUAL "c++" )
        #target_compile_definitions( feelpp PUBLIC _LIBCPP_ENABLE_CXX17_REMOVED_AUTO_PTR )
        target_compile_definitions( feelpp PUBLIC _LIBCPP_ENABLE_CXX17_REMOVED_FEATURES )
      endif()
    endif()
  endif()

if ( TARGET Boost::headers ) # from cmake >3.15
  target_link_libraries( feelpp PUBLIC Boost::headers )
elseif( TARGET Boost::boost )
  target_link_libraries( feelpp PUBLIC Boost::boost )
else()
  target_include_directories( feelpp PUBLIC ${Boost_INCLUDE_DIR} )
endif()
#target_link_libraries( feelpp PUBLIC ${Boost_LIBRARIES} )
foreach(COMPONENT ${BOOST_COMPONENTS_REQUIRED})
  target_link_libraries( feelpp PUBLIC Boost::${COMPONENT} )
endforeach()
endif()

if ( FFTW_FOUND )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_FFTW )
  target_include_directories( feelpp PRIVATE ${FFTW_INCLUDES} )
  target_link_libraries( feelpp PUBLIC ${FFTW_LIBRARIES} )
endif()
if ( ANN_FOUND )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_ANN FEELPP_HAS_ANN_H )
  target_include_directories( feelpp PUBLIC ${ANN_INCLUDE_DIR} )
  target_link_libraries( feelpp PUBLIC ${ANN_LIBRARIES} )
endif()
if ( GLPK_FOUND )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_GLPK FEELPP_HAS_GLPK_H )
  target_include_directories( feelpp PUBLIC ${GLPK_INCLUDE_DIR} )
  target_link_libraries( feelpp PUBLIC ${GLPK_LIBRARIES} )
endif()
if ( GOOGLE_PERFTOOLS_FOUND )
  target_include_directories( feelpp PRIVATE ${GOOGLE_PERFTOOLS_INCLUDE_DIR} )
  target_link_libraries( feelpp PUBLIC ${TCMALLOC_LIBRARIES} )
endif()
if (PETSC_FOUND)
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_PETSC FEELPP_HAS_PETSC_H )
  target_include_directories( feelpp PUBLIC ${PETSC_PACKAGE_INCLUDES} ${PETSC_INCLUDE_DIR} ${PETSC_INCLUDE_CONF})
  target_link_libraries( feelpp PUBLIC ${PETSC_LIBRARIES} ${PETSC_PACKAGE_LIBS} )
endif()
if (SLEPC_FOUND)
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_SLEPC )
  target_include_directories( feelpp PUBLIC ${SLEPC_INCLUDE_DIR} )
  target_link_libraries( feelpp PUBLIC ${SLEPC_LIBRARIES} )
endif()
if (GSL_FOUND)
  target_include_directories( feelpp PUBLIC ${GSL_INCLUDE_DIRS} )
  target_link_libraries( feelpp PUBLIC ${GSL_LIBRARIES} )
endif()

# Gmsh
if ( GMSH_FOUND )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_GMSH FEELPP_HAS_GMSH_H FEELPP_HAS_GMSH_LIBRARY GMSH_EXECUTABLE=${GMSH_EXECUTABLE} )
  if ( FEELPP_HAS_GMSH_API )
    target_compile_definitions( feelpp PUBLIC FEELPP_HAS_GMSH_API )
    if ( GMSH_LIBRARY_VERSION VERSION_LESS 4.2 )
      target_compile_definitions( feelpp PUBLIC GMSH_API_VERSION_MAJOR=${GMSH_MAJOR_VERSION} GMSH_API_VERSION_MINOR=${GMSH_MINOR_VERSION} )
    endif()
  endif()
  if ( FEELPP_HAS_MED )
    target_compile_definitions( feelpp PUBLIC FEELPP_HAS_MED FEELPP_HAS_GMSH_HAS_MED )
  endif()
  if ( FEELPP_HAS_GMSH_API )
    if ( NOT FEELPP_GMSH_API_MISSING_HEADERS )
      target_include_directories( feelpp PUBLIC ${GMSH_INCLUDE_DIR} )
    endif()
  else()
    target_include_directories( feelpp BEFORE PUBLIC ${GMSH_INCLUDE_DIR} )  # fix conflicts (on macos homebrew) with some Python sources
  endif()
  if ( GL2PS_LIBRARY )
    if ( GL_LIBRARY AND FEELPP_ENABLE_OPENGL )
      target_link_libraries( feelpp PUBLIC ${GL2PS_LIBRARY} ${GL_LIBRARY})
    else()
      target_link_libraries( feelpp PUBLIC ${GL2PS_LIBRARY})
    endif()
  endif()
  # GMSH_EXTERNAL_LIBRARIES may contain the med / cgns libs if they have been found
  target_link_libraries( feelpp PUBLIC ${GMSH_LIBRARIES} ${GMSH_EXTERNAL_LIBRARIES})
endif( GMSH_FOUND )

# Paraview
if ( ParaView_FOUND )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_PARAVIEW FEELPP_HAS_VTK )
  target_compile_definitions( feelpp PUBLIC VTK_HAS_PARALLEL=${VTK_HAS_PARALLEL} )
  if ( FEELPP_VTK_INSITU_ENABLED )
    target_compile_definitions( feelpp PUBLIC FEELPP_VTK_INSITU_ENABLED )
  endif()
  target_include_directories( feelpp PUBLIC ${ParaView_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})
  target_link_libraries( feelpp PUBLIC ${ParaView_LIBRARIES} ${VTK_LIBRARIES})
elseif ( VTK_FOUND )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_VTK )
  target_compile_definitions( feelpp PUBLIC VTK_HAS_PARALLEL=${VTK_HAS_PARALLEL} )
  target_include_directories( feelpp PUBLIC ${VTK_INCLUDE_DIRS})
  target_link_libraries( feelpp PUBLIC ${VTK_LIBRARIES})
endif()
if ( FEELPP_HAS_LIBXML2 )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_LIBXML2 )
  target_include_directories( feelpp PUBLIC ${LIBXML2_INCLUDE_DIR} )
  target_link_libraries( feelpp PUBLIC ${LIBXML2_LIBRARIES})
endif()
if ( FEELPP_HAS_ACUSIM )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_ACUSIM )
  target_include_directories( feelpp    PUBLIC ${ACUSIM_INCLUDE_DIR})
  target_link_libraries( feelpp PUBLIC ${ACUSIM_LIBRARIES} )
endif()
if (OMC_FOUND)
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_OMC )
#  target_include_directories( feelpp PUBLIC ${OMC_INCLUDE_DIR} )
#  target_link_libraries( feelpp PUBLIC ${OMC_LIBRARIES} )
endif()
if ( OPENTURNS_FOUND )
  target_compile_definitions( feelpp PUBLIC FEELPP_HAS_OPENTURNS )
  target_include_directories( feelpp PUBLIC ${OpenTURNS_INCLUDE_DIRS} )
  target_link_libraries( feelpp PUBLIC ${OpenTURNS_LIBRARIES} )
endif()

get_property(INC_DIR_LIST TARGET feelpp PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "[feelpp]: include_dir ${INC_DIR_LIST}")
get_property(INC_DIR_LIST TARGET feelpp PROPERTY PRIVATE_INCLUDE_DIRECTORIES)
message(STATUS "[feelpp]: private einclude_dir ${INC_DIR_LIST}")
get_property(INC_DIR_LIST TARGET feelpp PROPERTY PRIVATE_LINK_LIBRARIES)
message(STATUS "[feelpp]: private libs ${INC_DIR_LIST}")
get_property(INC_DIR_LIST TARGET feelpp PROPERTY PUBLIC_LINK_LIBRARIES)
message(STATUS "[feelpp]: public libs ${INC_DIR_LIST}")
get_property(INC_DIR_LIST TARGET feelpp PROPERTY LINK_LIBRARIES)
message(STATUS "[feelpp]: libs ${INC_DIR_LIST}")

#target_link_libraries(feelpp ${FEELPP_LIBRARIES})
if( FEELPP_ENABLE_PCH )
    add_precompiled_header( feelpp )
endif()
include(GNUInstallDirs)
# Install feelpp
install(TARGETS feelpp EXPORT FeelppTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
install(EXPORT FeelppTargets
  FILE FeelppTargets.cmake
  NAMESPACE Feelpp::
  DESTINATION share/feelpp/feel/cmake/modules
  )

#INSTALL(TARGETS feelpp DESTINATION lib/ COMPONENT Libs EXPORT feelpp-export-targets)
#install(EXPORT feelpp-export-targets DESTINATION share/feelpp/feel/cmake/modules)

# Handle static library if needed
# if(FEELPP_ENABLE_BUILD_STATIC)
#     add_library(feelpp-static STATIC options.cpp ${FEELPP_SRCS})
#     set_target_properties(feelpp-static PROPERTIES OUTPUT_NAME "feelpp")
#     set_target_properties(feelpp-static PROPERTIES VERSION ${FEELPP_SHARED_VERSION})
#     set_property(TARGET feelpp-static PROPERTY LABELS feel++)

#     target_link_libraries(feelpp-static ${FEELPP_LIBRARIES})

#     # Install feelpp
#     INSTALL(TARGETS feelpp-static DESTINATION lib/ COMPONENT Bin)
# endif()

include(feelpp.lib.install)
