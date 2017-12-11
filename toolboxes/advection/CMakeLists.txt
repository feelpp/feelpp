
feelpp_add_application( advection_2d SRCS main_adv.cpp
  LINK_LIBRARIES
  #feelpp_model_advection
  feelpp_model_advection2dP1G1DRP1
  feelpp_model_advection2dP2G1DRP2
  feelpp_model_advection2dP1VecG1DRP1
  feelpp_model_advection2dP2VecG1DRP2
  DEFS FEELPP_DIM=2
)

#feelpp_add_application( application_fluid_3d SRCS main_fm.cpp
  #LINK_LIBRARIES
  #feelpp_model_fluidmec3dP2P1G1
  #DEFS FEELPP_DIM=3
#)

feelpp_add_testcase(ring PREFIX feelpp_advection CATEGORY advection)
feelpp_add_testcase(square PREFIX feelpp_advection CATEGORY advection)

#feelpp_add_application( application_fm_2d SRCS main_fm.cpp LINK_LIBRARIES feelpp_model_fluidmechanics )

#install(TARGETS feelpp_application_fluid_2d feelpp_application_fluid_3d 
        #RUNTIME DESTINATION bin/ 
        #COMPONENT ModelApplications)
