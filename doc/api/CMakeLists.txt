
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Set to "ON" if you want the LaTeX documentation to be built.
set(DOXYFILE_LATEX 0)
include(UseDoxygen OPTIONAL)


configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/latexheader.tex.in
  ${CMAKE_CURRENT_BINARY_DIR}/latexheader.tex
)

