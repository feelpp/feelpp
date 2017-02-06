feelpp_add_application( thermoelectric_2d SRCS main_thermoelectric.cpp
  LINK_LIBRARIES
  feelpp_model_thermoelectric
  DEFS FEELPP_DIM=2
  INSTALL
  )

feelpp_add_application( thermoelectric_3d SRCS main_thermoelectric.cpp
  LINK_LIBRARIES
  feelpp_model_thermoelectric
  DEFS FEELPP_DIM=3
  INSTALL
  )

feelpp_add_testcase(test PREFIX thermoelectric CATEGORY ThermoElectric)
feelpp_add_testcase(ElectroMagnets PREFIX thermoelectric CATEGORY ThermoElectric)
