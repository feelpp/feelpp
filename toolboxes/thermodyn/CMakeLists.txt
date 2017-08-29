
feelpp_add_application( thermodyn_2d SRCS main_td.cpp
  LINK_LIBRARIES feelpp_model_thermodyn2dP1G1
  DEFS FEELPP_DIM=2
  INSTALL NO_TEST
  )

feelpp_add_application( thermodyn_3d SRCS main_td.cpp
  LINK_LIBRARIES feelpp_model_thermodyn3dP1G1
  DEFS FEELPP_DIM=3
  INSTALL NO_TEST
  )

feelpp_add_testcase(thermo2d CATEGORY HeatTransfer )
feelpp_add_testcase(Building CATEGORY HeatTransfer )
