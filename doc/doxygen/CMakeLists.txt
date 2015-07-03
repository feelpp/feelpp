#Include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )                                                 
#add_custom_target(doxygen)                                                                         
if(FEELPP_ENABLE_DOXYGEN)
find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it correctly")
  endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in  ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/footer.html  ${PROJECT_BINARY_DIR}/footer.html @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/header.html  ${PROJECT_BINARY_DIR}/header.html @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/layout.xml   ${PROJECT_BINARY_DIR}/layout.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/feel.css     ${PROJECT_BINARY_DIR}/feel.css @ONLY)

add_custom_command(
                  OUTPUT ${PROJECT_BINARY_DIR}/html/index.html
                         ${PROJECT_BINARY_DIR}/latex/Makefile
                         ${PROJECT_BINARY_DIR}/latex/refman.pdf

                  COMMAND  ${DOXYGEN}
                   ARGS     ${PROJECT_BINARY_DIR}/Doxyfile

                   DEPENDS  ${PROJECT_BINARY_DIR}/Doxyfile
                            ${PROJECT_BINARY_DIR}/html/index.html
                            ${PROJECT_BINARY_DIR}/latex/Makefile)

add_custom_target (doxy ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
            DEPENDS ${PROJECT_BINARY_DIR}/Doxyfile
                  # ${CMAKE_CURRENT_DIR}/html/index.html                                            
                 # ${CMAKE_CURRENT_BINARY_DIR}/latex/Makefile                                       
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM)

endif()
