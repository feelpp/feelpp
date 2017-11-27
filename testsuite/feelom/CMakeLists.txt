om_add_model( test_om SRCS test_om.mo )
feelpp_add_application( test_om SRCS test_om.cpp LINK_LIBRARIES ${OMLIB_DIR}/libtest_om.so )
#target_link_libraries( feelpp_test_test_om test_om.so )
