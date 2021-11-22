
feelpp_add_application( heat SRCS main_ht.cpp
  LINK_LIBRARIES
  # feelpp_toolbox_heat_lib_2dP1G1
  feelpp_toolbox_heat_lib
  INSTALL TESTS
  PROJECT toolbox
  MAN heat
  )

feelpp_add_testcase(cases PREFIX feelpp_toolbox_heat CATEGORY toolboxes/heat)

