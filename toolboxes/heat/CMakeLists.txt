
feelpp_add_application( heat SRCS main_ht.cpp
  LINK_LIBRARIES
  # feelpp_toolbox_heat_lib_2dP1G1
  feelpp_toolbox_heat_lib
  INSTALL TESTS
  )

feelpp_add_testcase(thermo2d CATEGORY Heat )
feelpp_add_testcase(Building CATEGORY Heat )
