
feelpp_add_application( heat_2d SRCS main_ht.cpp
  LINK_LIBRARIES
  # feelpp_toolbox_heat_lib_2dP1G1
  feelpp_toolbox_heat_lib
  DEFS FEELPP_DIM=2
  INSTALL
  )

feelpp_add_application( heat_3d SRCS main_ht.cpp
  LINK_LIBRARIES
  # feelpp_toolbox_heat_lib_3dP1G1
  feelpp_toolbox_heat_lib
  DEFS FEELPP_DIM=3
  INSTALL
  )

feelpp_add_testcase(thermo2d CATEGORY Heat )
feelpp_add_testcase(Building CATEGORY Heat )
