
# genLibThermoDynamics(
#   DIM 2
#   T_ORDER 1
#   GEO_ORDER 1
# )

feelpp_add_application( application_thermodyn_2d SRCS main_td.cpp LINK_LIBRARIES feelpp_model_thermodyn2dP1G1 )
#feelpp_add_application( application_thermodyn SRCS main_td.cpp LINK_LIBRARIES feelpp_model_thermodynamics )
