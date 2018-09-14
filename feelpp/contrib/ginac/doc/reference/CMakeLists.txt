set(VERSION ${GINAC_VERSION})
set(top_srcdir ${CMAKE_SOURCE_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DoxyfileHTML.in
	       ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileHTML
	       @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DoxyfilePDF.in
	       ${CMAKE_CURRENT_BINARY_DIR}/DoxyfilePDF
	       @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DoxyfileTEX.in
	       ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileTEX
	       @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfooter.in
	       ${CMAKE_CURRENT_BINARY_DIR}/Doxyfooter
	       @ONLY)

add_custom_target(html_dox
	COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileHTML
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "DOXYGEN DoxyfileHTML")
add_dependencies(html html_dox)
add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pdflatex/refman.tex
	COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/DoxyfilePDF
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "DOXYGEN DoxyfilePDF")
add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pdflatex/reference.tex
	COMMAND ${CMAKE_COMMAND} -E copy refman.tex reference.tex
	DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pdflatex/refman.tex
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pdflatex)


if (LATEX_FOUND)
	pdflatex_process(${CMAKE_CURRENT_BINARY_DIR}/pdflatex/reference.tex)
	add_custom_target(pdf_dox DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pdflatex/reference.pdf)
	add_dependencies(pdf pdf_dox)
endif()


