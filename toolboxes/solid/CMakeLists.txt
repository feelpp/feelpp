
feelpp_add_application( solid_2d SRCS main_sm.cpp
  LINK_LIBRARIES
  feelpp_toolbox_solid_lib_2dP1G1
  feelpp_toolbox_solid_lib_2dP2G1
  DEFS FEELPP_DIM=2
  INSTALL
 )

feelpp_add_application( solid_3d SRCS main_sm.cpp
  LINK_LIBRARIES
  feelpp_toolbox_solid_lib_3dP1G1
  feelpp_toolbox_solid_lib_3dP2G1
  DEFS FEELPP_DIM=3
  INSTALL
 )

feelpp_add_application( stress_2d SRCS stress.cpp
  LINK_LIBRARIES
  feelpp_toolbox_solid_lib_2dP1G1
  feelpp_toolbox_solid_lib_2dP2G1
  DEFS FEELPP_DIM=2
 )

feelpp_add_application( stress_3d SRCS stress.cpp
  LINK_LIBRARIES
  feelpp_toolbox_solid_lib_3dP1G1
  feelpp_toolbox_solid_lib_3dP2G1
  DEFS FEELPP_DIM=3
 )

feelpp_add_application( solenoid_3d SRCS stress.cpp
  LINK_LIBRARIES
  feelpp_toolbox_solid_lib_3dP1G1
  feelpp_toolbox_solid_lib_3dP2G1
  DEFS FEELPP_DIM=3 ELASTICITY_ANA
 )


#feelpp_add_application( application_solid SRCS main_sm.cpp LINK_LIBRARIES feelpp_model_solidmechanics )

feelpp_add_testcase(NAFEMS-LE1 CATEGORY CSM)
feelpp_add_testcase(NAFEMS-LE10 CATEGORY CSM)
feelpp_add_testcase(Solenoid CATEGORY CSM)
feelpp_add_testcase(TurekHron PREFIX solid CATEGORY CSM)
feelpp_add_testcase(torsionbar CATEGORY CSM)
feelpp_add_testcase(cantilever CATEGORY CSM)
