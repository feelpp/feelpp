find_package(Doxygen)
find_program(FIG2DEV fig2dev)
find_program(TEXI2DVI texi2dvi)
find_package(LATEX)
if (PDFLATEX_COMPILER AND MAKEINDEX_COMPILER)
	set(LATEX_FOUND true)
endif()

macro(add_info_file thename)
	set(${thename}_INFO ${CMAKE_CURRENT_BINARY_DIR}/${thename}.info)
	set(${thename}_HTML ${CMAKE_CURRENT_BINARY_DIR}/${thename}.html)
	set(${thename}_TEXINFO ${CMAKE_CURRENT_SOURCE_DIR}/${thename}.texi)
	add_custom_command(
		OUTPUT ${${thename}_INFO}
		COMMAND ${MAKEINFO} --no-split -o ${${thename}_INFO} ${${thename}_TEXINFO}
		DEPENDS ${${thename}_TEXINFO}
		COMMENT "MAKEFINO ${thename}.texi"
		VERBATIM)
	add_custom_command(
		OUTPUT ${${thename}_HTML}
		COMMAND ${MAKEINFO} --html --no-split -o ${${thename}_HTML} ${${thename}_TEXINFO}
		DEPENDS ${${thename}_TEXINFO} ${${thename}_HTML_EXTRA_DEPS}
		COMMENT "MAKEINFOHTML ${thename}.texi"
		VERBATIM)
	add_custom_target(${thename}_info ALL DEPENDS ${${thename}_INFO})
	add_custom_target(${thename}_html DEPENDS ${${thename}_HTML})
	add_dependencies(info ${thename}_info)
	add_dependencies(html ${thename}_html)
	install(FILES ${${thename}_INFO} DESTINATION "${SHARE_INSTALL_PREFIX}/info")
endmacro()

macro(pdflatex_process texfile)
	get_filename_component(_dirname "${texfile}" PATH)
	get_filename_component(_basename "${texfile}" NAME_WE)
	set(_idx ${_dirname}/${_basename}.idx)
	set(_ind ${_dirname}/${_basename}.ind)
	set(_pdf ${_dirname}/${_basename}.pdf)
	set(_fixupind ${CMAKE_SOURCE_DIR}/scripts/fixupind.py)
	add_custom_command(
		OUTPUT ${_idx}
		COMMAND ${PDFLATEX_COMPILER} ${texfile}
		WORKING_DIRECTORY ${_dirname}
		DEPENDS ${texfile}
		COMMENT "PDFLATEX ${_basename}.tex (1)")
	add_custom_command(
		OUTPUT ${_ind}
		COMMAND ${MAKEINDEX_COMPILER} ${_idx}
		COMMAND python ${_fixupind} ${_idx}
		WORKING_DIRECTORY ${_dirname}
		DEPENDS ${texfile} ${_idx}
		COMMENT "MAKEINDEX ${_basename}.idx")
	add_custom_command(
		OUTPUT ${_pdf}
		COMMAND ${PDFLATEX_COMPILER} ${texfile}
		WORKING_DIRECTORY ${_dirname}
		DEPENDS ${texfile} ${_ind}
		COMMENT "PDFLATEX ${_basename}.tex (2)")
endmacro()

add_subdirectory(tutorial)
add_subdirectory(examples)
if (DOXYGEN_FOUND)
	add_subdirectory(reference)
endif()

