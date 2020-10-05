feelpp_add_application( benchmark_fast_marching_2d SRCS main_fast_marching.cpp
  LINK_LIBRARIES
  feelpp_toolbox_levelset_lib_2dP1P1G1
  feelpp_toolbox_levelset_lib_2dP2P2G1
  DEFS FEELPP_DIM=2 LEVELSET_PN_ORDER=${FEELPP_MODELS_LEVELSET_PN_ORDER}
  INSTALL
)

feelpp_add_application( benchmark_fast_marching_3d SRCS main_fast_marching.cpp
  LINK_LIBRARIES
  feelpp_toolbox_levelset_lib_3dP1P1G1
  feelpp_toolbox_levelset_lib_3dP2P2G1
  DEFS FEELPP_DIM=3 LEVELSET_PN_ORDER=${FEELPP_MODELS_LEVELSET_PN_ORDER}
  INSTALL
)


