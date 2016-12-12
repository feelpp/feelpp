
feelpp_add_application( application_thermoelectric_2d SRCS main_thermoelectric.cpp
  LINK_LIBRARIES
  feelpp_model_thermoelectric2dP1G1
  feelpp_model_thermoelectric2dP2G1
  feelpp_model_thermoelectric2dP3G1
  DEFS FEELPP_DIM=2
  )

feelpp_add_application( application_thermoelectric_3d SRCS main_thermoelectric.cpp
  LINK_LIBRARIES
  feelpp_model_thermoelectric3dP1G1
  feelpp_model_thermoelectric3dP2G1
  feelpp_model_thermoelectric3dP3G1
  DEFS FEELPP_DIM=3
  )

feelpp_add_testcase(test PREFIX thermoelectric CATEGORY thermoelectric)
