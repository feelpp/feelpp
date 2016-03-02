
feelpp_add_application( application_thermodyn_2d SRCS main_td.cpp
  LINK_LIBRARIES feelpp_model_thermodyn2dP1G1
  DEFS FEELPP_DIM=2
  )

feelpp_add_application( application_thermodyn_3d SRCS main_td.cpp
  LINK_LIBRARIES feelpp_model_thermodyn3dP1G1
  DEFS FEELPP_DIM=3
  )
