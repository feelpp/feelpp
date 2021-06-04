
feelpp_add_application( advection_2d SRCS main_adv.cpp
  LINK_LIBRARIES
  #feelpp_model_advection
  feelpp_toolbox_advection_lib_2dP1G1DRP1
  feelpp_toolbox_advection_lib_2dP2G1DRP2
  feelpp_toolbox_advection_lib_2dP1VecG1DRP1
  feelpp_toolbox_advection_lib_2dP2VecG1DRP2
  DEFS FEELPP_DIM=2
  INSTALL
)

feelpp_add_testcase(cases PREFIX feelpp_toolbox_advection CATEGORY toolboxes/advection)

