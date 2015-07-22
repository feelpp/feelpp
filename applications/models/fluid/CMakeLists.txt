
feelpp_add_application( application_fluid_2d SRCS main_fm.cpp
  LINK_LIBRARIES
  feelpp_model_fluidmec2dP2P1G1
  feelpp_model_fluidmec2dP1P1G1
  #feelpp_model_fluidmec2dP2P1cG1DVP0c
  #feelpp_model_fluidmec2dP1P1cG1DVP0c
  DEFS FEELPP_DIM=2
)

#feelpp_add_application( application_fm_2d SRCS main_fm.cpp LINK_LIBRARIES feelpp_model_fluidmechanics )
