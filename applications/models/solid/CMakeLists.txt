
feelpp_add_application( application_solid_2d SRCS main_sm.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec2dP1G1
  feelpp_model_solidmec2dP2G1
  DEFS FEELPP_DIM=2
 )

feelpp_add_application( application_solid_3d SRCS main_sm.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec3dP1G1
  feelpp_model_solidmec3dP2G1
  DEFS FEELPP_DIM=3
 )

feelpp_add_application( application_stress_2d SRCS stress.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec2dP1G1
  feelpp_model_solidmec2dP2G1
  DEFS FEELPP_DIM=2
 )

feelpp_add_application( application_stress_3d SRCS stress.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec3dP1G1
  feelpp_model_solidmec3dP2G1
  DEFS FEELPP_DIM=3
 )

feelpp_add_application( application_solenoid_3d SRCS stress.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec3dP1G1
  feelpp_model_solidmec3dP2G1
  DEFS FEELPP_DIM=3 ELASTICITY_ANA
 )


#feelpp_add_application( application_solid SRCS main_sm.cpp LINK_LIBRARIES feelpp_model_solidmechanics )

install(TARGETS feelpp_application_solid_2d feelpp_application_solid_3d 
        RUNTIME DESTINATION bin/ 
        COMPONENT ModelApplications)

