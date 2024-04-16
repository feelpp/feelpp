find_program(MAKEINFO makeinfo)
find_program(TEXI2DVI texi2dvi)
find_package(LATEX)

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
	add_custom_target(${thename}_html ALL DEPENDS ${${thename}_HTML})
	add_dependencies(info ${thename}_info)
	add_dependencies(html ${thename}_html)
endmacro()

macro(add_doc_format src fmt)
	set(_out ${CMAKE_CURRENT_BINARY_DIR}/${src}.${fmt})
	set(_src ${CMAKE_CURRENT_SOURCE_DIR}/${src}.texi)
	add_custom_command(
		OUTPUT ${_out}
		COMMAND ${TEXI2DVI} --${fmt} --batch -o ${_out} ${_src}
		DEPENDS ${_src} ${_${src}_tutorial_figures_${fmt}}
		COMMENT "TEXI2DVI ${src}.texi => ${src}.${fmt}"
		VERBATIM)
	add_custom_target(${fmt}_${src}_tutorial DEPENDS ${_out})
	add_dependencies(${fmt} ${fmt}_${src}_tutorial)
endmacro()

if (MAKEINFO)
	add_info_file(cln)
	if (PDFLATEX_COMPILER)
		add_doc_format(cln pdf)
	endif()
endif()

