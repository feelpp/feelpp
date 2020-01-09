feelpp_add_application( maxwell SRCS main_maxwell.cpp
  LINK_LIBRARIES feelpp_toolbox_maxwell_lib
  INSTALL
  )

feelpp_add_application( maxwell_h1 SRCS main_maxwell.cpp
  LINK_LIBRARIES feelpp_toolbox_maxwell_h1_lib
  DEFS MAXWELL_H1
  INSTALL
  )

feelpp_add_application( cvg_maxwell SRCS cvg_maxwell.cpp #maxwell_test.cpp
  LINK_LIBRARIES feelpp_toolbox_maxwell_lib
  INSTALL
  )

feelpp_add_application( cvg_maxwell_h1 SRCS cvg_maxwell.cpp# maxwell_test.cpp
  LINK_LIBRARIES feelpp_toolbox_maxwell_h1_lib
  DEFS MAXWELL_H1
  INSTALL
  )

feelpp_add_testcase(test PREFIX maxwell CATEGORY maxwell )
feelpp_add_testcase(cvg PREFIX maxwell CATEGORY maxwell )
