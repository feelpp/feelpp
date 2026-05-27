
feelpp_add_application( electromagnetic SRCS main.cpp
  LINK_LIBRARIES
  feelpp_toolbox_electromagnetic_lib
  INSTALL #TESTS
 # MAN electromagnetic
  )

#feelpp_add_testcase(cases PREFIX feelpp_toolbox_electromagnetic CATEGORY toolboxes/electromagnetic)

