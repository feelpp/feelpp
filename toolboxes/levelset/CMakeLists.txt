
feelpp_add_application( levelset_2d SRCS main_ls.cpp
  LINK_LIBRARIES
  feelpp_model_levelset2dP1G1
  feelpp_model_levelset2dP2G1
  DEFS FEELPP_DIM=2 LEVELSET_PN_ORDER=${FEELPP_MODELS_LEVELSET_PN_ORDER}
)

feelpp_add_testcase(square PREFIX feelpp_levelset CATEGORY levelset )
feelpp_add_testcase(tube PREFIX feelpp_levelset CATEGORY levelset )

