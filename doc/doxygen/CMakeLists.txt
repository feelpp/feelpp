INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
if(FEELPP_ENABLE_DOXYGEN)
find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it correctly")
  endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in 
                       ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/footer.html ${PROJECT_BINARY_DIR}/footer.html @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/header.html ${PROJECT_BINARY_DIR}/header.html @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/layout.xml ${PROJECT_BINARY_DIR}/layout.xml @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/feel.css ${PROJECT_BINARY_DIR}/feel.css @ONLY)
   
add_custom_command(	
                  OUTPUT ${PROJECT_BINARY_DIR}/html/index.html
                         ${PROJECT_BINARY_DIR}/latex/Makefile 
                         ${PROJECT_BINARY_DIR}/latex/refman.tex)           
add_custom_target (doc ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
            SOURCES ${PROJECT_SOURCE_DIR}/doc/doxygen/Doxygen.in
            COMMENT "Generating API documentation with Doxygen" VERBATIM))  

  add_custom_target(html
           DEPENDS ${PROJECT_BINARY_DIR}/html/index.html
                    ${PROJECT_BINARY_DIR}/latex/Makefile)
endif()
