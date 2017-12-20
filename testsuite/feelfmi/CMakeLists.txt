if ( FEELPP_HAS_OMC AND FEELPP_HAS_FMILIB )
  add_ommodel( fmi2cs CLASS BouncingBall SRCS BouncingBall.mo TYPE cs VERS 2.0 )
  feelpp_add_test( fmi DEPS feelpp_add_ommodel_fmi2cs CFG test_fmi.cfg )
endif()
