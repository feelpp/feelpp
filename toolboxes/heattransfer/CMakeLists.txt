
feelpp_add_application( heattransfer_2d SRCS main_ht.cpp
  LINK_LIBRARIES feelpp_toolbox_heattransfer_lib_2dP1G1
  DEFS FEELPP_DIM=2
  INSTALL
  )

feelpp_add_application( heattransfer_3d SRCS main_ht.cpp
  LINK_LIBRARIES feelpp_toolbox_heattransfer_lib_3dP1G1
  DEFS FEELPP_DIM=3
  INSTALL
  )

feelpp_add_testcase(thermo2d CATEGORY HeatTransfer )
feelpp_add_testcase(Building CATEGORY HeatTransfer )
