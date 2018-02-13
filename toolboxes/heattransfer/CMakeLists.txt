
feelpp_add_application( heat_transfer_2d SRCS main_ht.cpp
  LINK_LIBRARIES feelpp_model_heattransfer2dP1G1
  DEFS FEELPP_DIM=2
  INSTALL
  )

feelpp_add_application( heat_transfer_3d SRCS main_ht.cpp
  LINK_LIBRARIES feelpp_model_heattransfer3dP1G1
  DEFS FEELPP_DIM=3
  INSTALL
  )

feelpp_add_testcase(thermo2d CATEGORY HeatTransfer )
feelpp_add_testcase(Building CATEGORY HeatTransfer )
