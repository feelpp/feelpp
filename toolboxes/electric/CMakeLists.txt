
feelpp_add_application( electric SRCS main_electric.cpp
  LINK_LIBRARIES
  # feelpp_toolbox_heat_lib_2dP1G1
  feelpp_toolbox_electric_lib
  INSTALL
  MAN electric
  )

feelpp_add_testcase(cases PREFIX feelpp_toolbox_electric CATEGORY toolboxes/electric)

