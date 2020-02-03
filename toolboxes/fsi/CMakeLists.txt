
feelpp_add_application( fsi_2d SRCS main_fsi.cpp
  LINK_LIBRARIES
  feelpp_toolbox_fsi_lib_2dP2P1G1_2dP1G1
  feelpp_toolbox_fsi_lib_2dP2P1G1_2dP2G1
  #feelpp_toolbox_fsi_lib_2dP1P1G1_2dP1G1
  DEFS FEELPP_DIM=2 FEELPP_GEO_ORDER=1
  INSTALL
  MAN fsi
  )

feelpp_add_application( fsi_3d SRCS main_fsi.cpp
  LINK_LIBRARIES
  feelpp_toolbox_fsi_lib_3dP2P1G1_3dP1G1
  feelpp_toolbox_fsi_lib_3dP2P1G1_3dP2G1
  DEFS FEELPP_DIM=3 FEELPP_GEO_ORDER=1
  INSTALL
  )

feelpp_add_application( fsi_2d_g2 SRCS main_fsi.cpp
  LINK_LIBRARIES
  feelpp_toolbox_fsi_lib_2dP2P1G2_2dP2G2
  DEFS FEELPP_DIM=2 FEELPP_GEO_ORDER=2
  INSTALL
  )

feelpp_add_application( fsi_3d_g2 SRCS main_fsi.cpp
  LINK_LIBRARIES
  feelpp_toolbox_fsi_lib_3dP2P1G2_3dP2G2
  DEFS FEELPP_DIM=3 FEELPP_GEO_ORDER=2
  INSTALL
  )


feelpp_add_testcase( TurekHron PREFIX feelmodels_fsi CATEGORY FSI)
feelpp_add_testcase( wavepressure2d PREFIX feelmodels_fsi CATEGORY FSI)      
feelpp_add_testcase( wavepressure3d PREFIX feelmodels_fsi CATEGORY FSI)      
feelpp_add_testcase( movingboundary PREFIX feelmodels_fsi CATEGORY FSI)
