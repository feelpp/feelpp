
feelpp_add_application( solid SRCS main_sm.cpp
  LINK_LIBRARIES
  feelpp_toolbox_solid_lib_2dP1G1
  feelpp_toolbox_solid_lib_2dP2G1
  feelpp_toolbox_solid_lib_3dP1G1
  feelpp_toolbox_solid_lib_3dP2G1
  INSTALL TESTS
  MAN solid
 )
if ( 0 )
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
endif()

feelpp_add_testcase(cases PREFIX feelpp_toolbox_solid CATEGORY toolboxes/solid)

add_dependencies(toolboxes feelpp_toolbox_solid)