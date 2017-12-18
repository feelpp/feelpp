set( SRC
  fmu.cpp
  fmumodel2.cpp
  )
file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
feel_append_src(feelfmi "${SRC};${HDR}")
