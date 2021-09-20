
feelpp_add_application( multifluid_2d SRCS main_mf.cpp
  LINK_LIBRARIES
  #feelpp_model_multifluid
  feelpp_toolbox_multifluid_lib_2dP2P1G1DVP1c_2dP1P2G1
  feelpp_toolbox_multifluid_lib_2dP2P1G1DVP1c_2dP2P2G1
  #feelpp_toolbox_multifluid_lib_2dP2P1G1DVP1c_2dP3P2G1
  DEFS FEELPP_DIM=2 FEELPP_GEO_ORDER=1 LEVELSET_PN_ORDER=${FEELPP_MODELS_LEVELSET_PN_ORDER}
  INSTALL
)

feelpp_add_application( multifluid_3d SRCS main_mf.cpp
  LINK_LIBRARIES
  #feelpp_model_multifluid
  feelpp_toolbox_multifluid_lib_3dP2P1G1DVP1c_3dP1P2G1
  feelpp_toolbox_multifluid_lib_3dP2P1G1DVP1c_3dP2P2G1
  #feelpp_toolbox_multifluid_lib_3dP2P1G1DVP1c_3dP3P2G1
  DEFS FEELPP_DIM=3 FEELPP_GEO_ORDER=1 LEVELSET_PN_ORDER=${FEELPP_MODELS_LEVELSET_PN_ORDER}
  INSTALL
)


feelpp_add_testcase(cases PREFIX feelpp_toolbox_multifluid CATEGORY toolboxes/multifluid)

