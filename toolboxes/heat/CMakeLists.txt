
feelpp_add_application( heat SRCS main_ht.cpp
  LINK_LIBRARIES
  # feelpp_toolbox_heat_lib_2dP1G1
  feelpp_toolbox_heat_lib
  INSTALL TESTS
  PROJECT toolbox
  MAN heat
  )

foreach( _testcase IN ITEMS "thermo2d" "Building" "h-convergence")
  feelpp_add_testcase(${_testcase} PREFIX feelpp_toolbox_heat CATEGORY toolboxes/heat)
endforeach()
