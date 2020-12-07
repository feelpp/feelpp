
feelpp_add_application( coefficientformpdes SRCS main.cpp
  LINK_LIBRARIES
  feelpp_toolbox_coefficientformpdes_lib
  #feelpp_toolbox_coefficientformpdes_2dG1
  INSTALL TESTS
  PROJECT toolbox
#  MAN coefficientformpdes
  )

# feelpp_add_application( coefficientformpdes_predefined SRCS main_predefined.cpp
#   LINK_LIBRARIES
#   feelpp_toolbox_coefficientformpdes_lib
#   #feelpp_toolbox_coefficientformpdes_2dG1
#   INSTALL #TESTS
#   PROJECT toolbox
# #  MAN coefficientformpdes
#   )
