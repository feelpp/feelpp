cmake_minimum_required(VERSION 2.8)
project(GKlib)


set(GKLIB_PATH ${FEELPP_SOURCE_DIR}/feelpp/contrib/metis/GKlib)
message(STATUS "GKlib path ${GKLIB_PATH}")
include(GKlibSystem.cmake)

include_directories(".")
add_library(GKlib STATIC ${GKlib_sources})
if(UNIX)
  target_link_libraries(GKlib m)
endif(UNIX)

include_directories("test")
add_subdirectory("test")

install(TARGETS GKlib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${GKlib_includes} DESTINATION include)
