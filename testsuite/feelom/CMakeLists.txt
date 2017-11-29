om_add_model( test_om SRCS test_om.mo LIB TEST_MO_LIB )
feelpp_add_application( test_om SRCS test_om.cpp LINK_LIBRARIES ${TEST_MO_LIB} )
#target_link_libraries( feelpp_test_test_om test_om.so )
