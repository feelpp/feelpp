# Add this directory for internal users.

# Find sources.
file(GLOB metis_sources *.c)

foreach(file ${metis_sources})
  get_filename_component( FILE_EXT ${file} EXT )
  #message(STATUS "extension ${file} : ${FILE_EXT}")
#  if ( ${FILE_EXT} STREQUAL ".c" )
#    message(STATUS "set ${file} as c++")
#    set_source_files_properties(${file} PROPERTIES LANGUAGE CXX)
#  endif()
endforeach()

# Build libmetis.
add_library(feelpp_metis ${METIS_LIBRARY_TYPE} ${GKlib_sources} ${metis_sources})
target_include_directories(feelpp_metis PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
if(UNIX)
 target_link_libraries(feelpp_metis m)
endif()
if ("${METIS_LIBRARY_TYPE}" STREQUAL "SHARED" )
  set_property(TARGET feelpp_metis PROPERTY MACOSX_RPATH ON)
endif()


if(METIS_INSTALL)
  install(TARGETS feelpp_metis
    EXPORT feelpp-contrib-export-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
