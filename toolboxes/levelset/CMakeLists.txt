
feelpp_add_application( levelset_2d SRCS main_ls.cpp
  LINK_LIBRARIES
  feelpp_model_levelset2dP1G1
  feelpp_model_levelset2dP2G1
  #feelpp_model_fluidmec2dP1P1G1
  #feelpp_model_fluidmec2dP2P1cG1DVP0c
  #feelpp_model_fluidmec2dP1P1cG1DVP0c
  DEFS FEELPP_DIM=2
)

feelpp_add_testcase(square PREFIX feelpp_levelset CATEGORY levelset )
feelpp_add_testcase(tube PREFIX feelpp_levelset CATEGORY levelset )

#feelpp_add_application( application_fm_2d SRCS main_fm.cpp LINK_LIBRARIES feelpp_model_fluidmechanics )

#install(TARGETS feelpp_application_fluid_2d feelpp_application_fluid_3d 
        #RUNTIME DESTINATION bin/ 
        #COMPONENT ModelApplications)
