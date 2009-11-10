IF ( ENABLE_INSTANTIATION_MODE )
  SET( LIFEDISCR_INST_SRCS mesh1.cpp mesh2.cpp mesh3.cpp partitioner.cpp ale.cpp meshhighorder.cpp )
endif()

add_library( lifediscr
  bdf2.cpp oseendata.cpp ${LIFEDISCR_INST_SRCS}
  )

target_link_libraries (lifediscr lifemesh lifealg)
