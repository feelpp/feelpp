
feelpp_add_application( coefficientformpdes SRCS main.cpp
  LINK_LIBRARIES
  feelpp_toolbox_coefficientformpdes_lib
  #feelpp_toolbox_coefficientformpdes_2dG1
  #feelpp_toolbox_coefficientformpdes_3dG1
  INSTALL TESTS
  PROJECT toolbox
#  MAN coefficientformpdes
  )

foreach( _testcase IN ITEMS "adr" "bratu" "heat" "mesh-motion" "predator-prey" "thermoelectric" "cahn-hilliard")
  feelpp_add_testcase(${_testcase} PREFIX feelpp_toolbox_coefficientformpdes CATEGORY toolboxes/coefficientformpdes)
endforeach()
