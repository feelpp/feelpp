
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

include(UseDoxygen)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/latexheader.tex.in
  ${CMAKE_CURRENT_BINARY_DIR}/latexheader.tex
)


