get_filename_component(FeelppContrib_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(CMakeFindDependencyMacro)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR})
foreach( dep feelpp_gflags glog nlopt )
  set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${FEELPP_DIR}/share/feelpp/${dep}/cmake )
endforeach()

find_dependency( feelpp_gflags REQUIRED )
find_dependency( glog REQUIRED )
find_dependency( Eigen3 REQUIRED )
find_dependency( pybind11 )
find_dependency( NLopt )

if(NOT TARGET Feelpp::feelpp_contrib)
  include("${FeelppContrib_CMAKE_DIR}/feelpp-contrib-export-targets.cmake")
endif()
