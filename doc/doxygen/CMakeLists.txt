if(FEELPP_ENABLE_DOXYGEN)
find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it correctly")
  endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in  ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)

add_custom_command(
                  OUTPUT ${PROJECT_BINARY_DIR}/html/index.html
                        
                  COMMAND  ${DOXYGEN}
                   ARGS     ${PROJECT_BINARY_DIR}/Doxyfile
                   DEPENDS  ${PROJECT_BINARY_DIR}/Doxyfile
add_custom_target (doxy ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
            DEPENDS ${PROJECT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM)

endif()
