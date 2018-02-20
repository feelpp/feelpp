feelpp_add_application( thermoelectric_2d SRCS main_thermoelectric.cpp
  LINK_LIBRARIES
  feelpp_toolbox_thermoelectric_lib
  #feelpp_toolbox_thermoelectric_lib_2dP1G1_2dP1G1
  DEFS FEELPP_DIM=2
  INSTALL
  )

feelpp_add_application( thermoelectric_3d SRCS main_thermoelectric.cpp
  LINK_LIBRARIES
  feelpp_toolbox_thermoelectric_lib
  #feelpp_toolbox_thermoelectric_lib_3dP1G1_3dP1G1
  DEFS FEELPP_DIM=3
  INSTALL
  )

feelpp_add_testcase(test PREFIX thermoelectric CATEGORY ThermoElectric)
feelpp_add_testcase(ElectroMagnets PREFIX thermoelectric CATEGORY ThermoElectric)
