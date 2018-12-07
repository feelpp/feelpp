configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.texi.cmake
	       ${CMAKE_CURRENT_BINARY_DIR}/version.texi @ONLY)

set(_ginac_tutorial_figures classhierarchy repnaive reppair repreal)
set(_ginac_tutorial_figures_pdf)
set(_ginac_tutorial_figures_png)
set(_ginac_tutorial_figures_eps)

macro(add_xfig_figure thename)
	set(${thename}_XFIG ${CMAKE_CURRENT_SOURCE_DIR}/${thename}.fig)
	set(${thename}_PDF ${CMAKE_CURRENT_BINARY_DIR}/${thename}.pdf)
	set(${thename}_EPS ${CMAKE_CURRENT_BINARY_DIR}/${thename}.eps)
	set(${thename}_PNG ${CMAKE_CURRENT_BINARY_DIR}/${thename}.png)
	list(APPEND _ginac_tutorial_figures_pdf ${${thename}_PDF})
	list(APPEND _ginac_tutorial_figures_png ${${thename}_PNG})
	list(APPEND _ginac_tutorial_figures_png ${${thename}_EPS})
	add_custom_command(
		OUTPUT ${${thename}_PDF}
		COMMAND ${FIG2DEV} -L pdf -m 0.9 ${${thename}_XFIG} ${${thename}_PDF}
		DEPENDS ${${thename}_XFIG}
		COMMENT "FIG2DEV ${thename}.xfig ==> ${thename}.pdf"
		VERBATIM)
	add_custom_command(
		OUTPUT ${${thename}_EPS}
		COMMAND ${FIG2DEV} -L eps -m 0.9 ${${thename}_XFIG} ${${thename}_EPS}
		DEPENDS ${${thename}_XFIG}
		COMMENT "FIG2DEV ${thename}.xfig ==> ${thename}.eps"
		VERBATIM)
	add_custom_command(
		OUTPUT ${${thename}_PNG}
		COMMAND ${FIG2DEV} -L png -m 0.9 ${${thename}_XFIG} ${${thename}_PNG}
		DEPENDS ${${thename}_XFIG}
		COMMENT "FIG2DEV ${thename}.xfig ==> ${thename}.png"
		VERBATIM)
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

if (TEXI2DVI AND FIG2DEV)
	foreach(_f ${_ginac_tutorial_figures})
		add_xfig_figure(${_f})
	endforeach()
	set(ginac_HTML_EXTRA_DEPS ${_ginac_tutorial_figures_png})
	add_doc_format(ginac pdf)
endif()

add_info_file(ginac)

