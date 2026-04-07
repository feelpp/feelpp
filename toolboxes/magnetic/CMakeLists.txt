
feelpp_add_application( magnetic SRCS main.cpp
  LINK_LIBRARIES
  feelpp_toolbox_magnetic_lib
  INSTALL
  MAN magnetic
  )

#feelpp_add_testcase(cases PREFIX feelpp_toolbox_magnetic CATEGORY toolboxes/magnetic)

