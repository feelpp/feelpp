
feelpp_add_application( solid_2d SRCS main_sm.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec2dP1G1
  feelpp_model_solidmec2dP2G1
  DEFS FEELPP_DIM=2
  INSTALL NO_TEST
 )

feelpp_add_application( solid_3d SRCS main_sm.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec3dP1G1
  feelpp_model_solidmec3dP2G1
  DEFS FEELPP_DIM=3
  INSTALL NO_TEST
 )

feelpp_add_application( stress_2d SRCS stress.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec2dP1G1
  feelpp_model_solidmec2dP2G1
  DEFS FEELPP_DIM=2 NO_TEST
 )

feelpp_add_application( stress_3d SRCS stress.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec3dP1G1
  feelpp_model_solidmec3dP2G1
  DEFS FEELPP_DIM=3 NO_TEST
 )

feelpp_add_application( solenoid_3d SRCS stress.cpp
  LINK_LIBRARIES
  feelpp_model_solidmec3dP1G1
  feelpp_model_solidmec3dP2G1
  DEFS FEELPP_DIM=3 ELASTICITY_ANA NO_TEST
 )


#feelpp_add_application( application_solid SRCS main_sm.cpp LINK_LIBRARIES feelpp_model_solidmechanics )

feelpp_add_testcase(NAFEMS-LE1 CATEGORY CSM)
feelpp_add_testcase(NAFEMS-LE10 CATEGORY CSM)
feelpp_add_testcase(Solenoid CATEGORY CSM)
feelpp_add_testcase(TurekHron PREFIX solid CATEGORY CSM)
feelpp_add_testcase(torsionbar CATEGORY CSM)
