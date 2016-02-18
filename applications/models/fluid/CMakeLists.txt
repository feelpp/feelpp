
feelpp_add_application( application_fluid_2d SRCS main_fm.cpp
  LINK_LIBRARIES
  feelpp_model_fluidmec2dP2P1G1
  #feelpp_model_fluidmec2dP1P1G1
  #feelpp_model_fluidmec2dP2P1cG1DVP0c
  #feelpp_model_fluidmec2dP1P1cG1DVP0c
  DEFS FEELPP_DIM=2
)

feelpp_add_application( application_fluid_3d SRCS main_fm.cpp
  LINK_LIBRARIES
  feelpp_model_fluidmec3dP2P1G1
  DEFS FEELPP_DIM=3
)

feelpp_add_testcase(cavity)
feelpp_add_testcase(TurekHron)
feelpp_add_testcase(TurekHronExt3d)

#feelpp_add_application( application_fm_2d SRCS main_fm.cpp LINK_LIBRARIES feelpp_model_fluidmechanics )

install(TARGETS feelpp_application_fluid_2d feelpp_application_fluid_3d 
        RUNTIME DESTINATION bin/ 
        COMPONENT ModelApplications)
