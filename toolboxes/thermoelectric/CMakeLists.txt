feelpp_add_application( thermoelectric SRCS main_thermoelectric.cpp
  LINK_LIBRARIES
  feelpp_toolbox_thermoelectric_lib
  #feelpp_toolbox_thermoelectric_lib_3dP1G1_3dP1G1
  #feelpp_toolbox_thermoelectric_lib_2dP1G1_2dP1G1
  INSTALL TESTS
  MAN thermoelectric
  )
if ( 0 )
feelpp_add_application( cvg_thermoelectric_2d SRCS cvg_thermoelectric.cpp
  LINK_LIBRARIES
  feelpp_toolbox_thermoelectric_lib
  #feelpp_toolbox_thermoelectric_lib_2dP1G1_2dP1G1
  DEFS FEELPP_DIM=2
  INSTALL
  )

feelpp_add_application( cvg_thermoelectric_3d SRCS cvg_thermoelectric.cpp
  LINK_LIBRARIES
  feelpp_toolbox_thermoelectric_lib
  #feelpp_toolbox_thermoelectric_lib_3dP1G1_3dP1G1
  DEFS FEELPP_DIM=3
  INSTALL
  )
endif()

feelpp_add_testcase(cases PREFIX feelpp_toolbox_thermoelectric CATEGORY toolboxes/thermoelectric)
