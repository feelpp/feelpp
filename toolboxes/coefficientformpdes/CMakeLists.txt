
feelpp_add_application( coefficientformpdes SRCS main.cpp
  LINK_LIBRARIES
  feelpp_toolbox_coefficientformpdes_lib
  #feelpp_toolbox_coefficientformpdes_2dG1
  #feelpp_toolbox_coefficientformpdes_3dG1
  INSTALL TESTS
  PROJECT toolbox
#  MAN coefficientformpdes
  )

feelpp_add_testcase(cases PREFIX feelpp_toolbox_cfpdes CATEGORY toolboxes/coefficientformpdes)

