
feelpp_add_application( levelset_2d SRCS main_ls.cpp
  LINK_LIBRARIES
  feelpp_toolbox_levelset_lib_2dP1G1
  feelpp_toolbox_levelset_lib_2dP2G1
  DEFS FEELPP_DIM=2
  INSTALL
)

feelpp_add_testcase(square PREFIX feelpp_levelset CATEGORY levelset )
feelpp_add_testcase(tube PREFIX feelpp_levelset CATEGORY levelset )
