
# genLibFSI(
#   DIM 2
#   FLUID_U_ORDER 2
#   FLUID_P_ORDER 1
#   #FLUID_P_CONTINUITY Discontinuous
#   FLUID_GEO_ORDER 1
#   SOLID_DISP_ORDER 1
#   SOLID_GEO_ORDER 1
# )

feelpp_add_application( application_fsi_2d SRCS main_fsi.cpp LINK_LIBRARIES feelpp_model_fsi_2dP2P1cG1DVP0c_2dP1G1DCLCst DEFS FEELPP_DIM=2)

feelpp_add_application( application_fsi_3d SRCS main_fsi.cpp LINK_LIBRARIES feelpp_model_fsi_3dP2P1cG1DVP0c_3dP1G1DCLCst DEFS FEELPP_DIM=3)

