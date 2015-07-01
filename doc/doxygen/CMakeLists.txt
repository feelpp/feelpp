INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR} )
if(FEELPP_ENABLE_DOXYGEN)
find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it correctly")
  endif()
  configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile.in 
                       ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/footer.html ${PROJECT_BINARY_DIR}/doc/doxygen/footer.html @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/header.html ${PROJECT_BINARY_DIR}/doc/doxygen/header.html @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/layout.xml ${PROJECT_BINARY_DIR}/doc/doxygen/layout.xml @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/feel.css ${PROJECT_BINARY_DIR}/doc/doxygen/feel.css @ONLY)
  
   
  add_custom_command(
<<<<<<< HEAD
	OUTPUT ${PROJECT_BINARY_DIR}/doc/html/index.html)
=======
		OUTPUT ${PROJECT_BINARY_DIR}/doc/html/index.html)
>>>>>>> 323057c6f3bb63026cea1fe29c141f2d10ad97b6
  add_custom_target (doc ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile
            SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxygen.in
            COMMENT "Generating API documentation with Doxygen" VERBATIM))  
            
  
<<<<<<< HEAD
  add_custom_target(html
	         DEPENDS ${PROJECT_BINARY_DIR}/doc/html/index.html	)
		           
=======
	add_custom_target(html
		         DEPENDS ${PROJECT_BINARY_DIR}/doc/html/index.html	)
	          
>>>>>>> 323057c6f3bb63026cea1fe29c141f2d10ad97b6
endif()
