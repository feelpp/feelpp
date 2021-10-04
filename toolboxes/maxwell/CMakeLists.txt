feelpp_add_application( maxwell_2d SRCS main_maxwell.cpp
  LINK_LIBRARIES feelpp_toolbox_maxwell_lib
  DEFS FEELPP_DIM=2
  INSTALL
  )

feelpp_add_application( maxwell_3d SRCS main_maxwell.cpp
  LINK_LIBRARIES feelpp_toolbox_maxwell_lib
  DEFS FEELPP_DIM=3
  INSTALL
  )

feelpp_add_application( cvg_maxwell_2d SRCS cvg_maxwell.cpp
  LINK_LIBRARIES feelpp_toolbox_maxwell_lib
  DEFS FEELPP_DIM=2
  INSTALL
  )

feelpp_add_application( cvg_maxwell_3d SRCS cvg_maxwell.cpp
  LINK_LIBRARIES feelpp_toolbox_maxwell_lib
  DEFS FEELPP_DIM=3
  INSTALL
  )

feelpp_add_application( biotsavart SRCS biotsavart_main.cpp INSTALL )
target_include_directories(feelpp_toolbox_biotsavart PRIVATE ${CMAKE_SOURCE_DIR})

feelpp_add_testcase(cases PREFIX feelpp_toolbox_maxwell CATEGORY toolboxes/maxwell)
