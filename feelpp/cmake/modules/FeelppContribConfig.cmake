get_filename_component(FeelppContrib_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(CMakeFindDependencyMacro)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR})
foreach( dep nlopt eigen3 )
  if ( EXISTS ${FEELPP_DIR}/share/feelpp/${dep}/cmake )
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${FEELPP_DIR}/share/feelpp/${dep}/cmake )
  endif()
endforeach()

find_dependency(Boost 1.65.0 REQUIRED COMPONENTS filesystem serialization mpi OPTIONAL_COMPONENTS system)
if (NOT TARGET Boost::system)
  add_library(Boost::system INTERFACE IMPORTED)
  set_property(TARGET Boost::system PROPERTY
    INTERFACE_COMPILE_DEFINITIONS "BOOST_ERROR_CODE_HEADER_ONLY;BOOST_SYSTEM_NO_LIB")
endif()

find_package(Eigen3 CONFIG QUIET)
if(NOT TARGET Eigen3::Eigen)
  find_dependency(Eigen3 REQUIRED)
endif()
if(NOT TARGET eigen AND TARGET Eigen3::Eigen)
  add_library(eigen INTERFACE IMPORTED)
  set_property(TARGET eigen PROPERTY INTERFACE_LINK_LIBRARIES Eigen3::Eigen)
endif()
if(NOT TARGET Eigen3::Eigen AND TARGET eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
  set_property(TARGET Eigen3::Eigen PROPERTY INTERFACE_LINK_LIBRARIES eigen)
endif()
find_dependency( pybind11 )
find_dependency( tabulate )
find_dependency( indicators )

# cln
find_package(PkgConfig REQUIRED)
pkg_search_module(CLN REQUIRED IMPORTED_TARGET "cln>=1.3.6")
message(STATUS "[feelpp] External CLN Includes: ${CLN_INCLUDE_DIRS}")
message(STATUS "[feelpp] External CLN Libraries: ${CLN_LIBRARIES}, ${CLN_LINK_LIBRARIES}")

if (CLN_FOUND AND NOT TARGET cln::cln)
  add_library(cln::cln INTERFACE IMPORTED)
  # Either forward to the pkg-config imported target:
  set_property(TARGET cln::cln PROPERTY INTERFACE_LINK_LIBRARIES PkgConfig::CLN)
  # And (optional) expose include dirs explicitly for IDEs:
  set_property(TARGET cln::cln PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CLN_INCLUDE_DIRS}")
endif()

if ( FEELPP_HAS_MMG )
  find_dependency( mmg )
endif()
if( FEELPP_HAS_PARMMG )
  find_dependency( parmmg )
endif()
find_dependency( range-v3 )
# find_dependency( specx )
find_dependency( CURL )
if ( EXISTS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules/swsConfig.cmake )
  find_dependency( sws )
endif()

find_dependency( NLopt )

# if ( FEELPP_HAS_SPECX )
#   find_dependency( specx )
# endif()
# 
# if ( FEELPP_HAS_EIGENRAND )
#   find_dependency( eigenrand )
# endif()


if(NOT TARGET Feelpp::feelpp_contrib)
  include("${FeelppContrib_CMAKE_DIR}/feelpp-contrib-export-targets.cmake")
endif()
