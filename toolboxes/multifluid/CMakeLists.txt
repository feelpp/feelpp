
feelpp_add_application( multifluid_2d SRCS main_mf.cpp
  LINK_LIBRARIES
  #feelpp_model_multifluid
  feelpp_toolbox_multifluid_lib_2dP2P1G1_2dP1G1
  feelpp_toolbox_multifluid_lib_2dP2P1G1_2dP2G1
  feelpp_toolbox_multifluid_lib_2dP2P1G1_2dP3G1
  DEFS FEELPP_DIM=2 FEELPP_GEO_ORDER=1
  INSTALL
)

feelpp_add_application( multifluid_3d SRCS main_mf.cpp
  LINK_LIBRARIES
  #feelpp_model_multifluid
  feelpp_toolbox_multifluid_lib_3dP2P1G1_3dP1G1
  feelpp_toolbox_multifluid_lib_3dP2P1G1_3dP2G1
  feelpp_toolbox_multifluid_lib_3dP2P1G1_3dP3G1
  DEFS FEELPP_DIM=3 FEELPP_GEO_ORDER=1
  INSTALL
)

feelpp_add_testcase(square PREFIX feelpp_multifluid CATEGORY multifluid)
feelpp_add_testcase(drop PREFIX feelpp_multifluid CATEGORY multifluid)
feelpp_add_testcase(vesicle_rest PREFIX feelpp_multifluid CATEGORY multifluid)
feelpp_add_testcase(vesicle_shear PREFIX feelpp_multifluid CATEGORY multifluid)
feelpp_add_testcase(sheared_sphere PREFIX feelpp_multifluid CATEGORY multifluid)
