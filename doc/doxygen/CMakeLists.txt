INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
if(FEELPP_ENABLE_DOXYGEN)
find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it correctly")
  endif()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in 
                       ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc ALL
             COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
                  COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif()
