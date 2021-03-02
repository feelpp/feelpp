
crb_add_library( toolbox-mor SRCS toolboxmor.cpp EXEC CRB_LIB_TOOLBOXMOR )

set(FEELPP_CRB_toolboxmor_LIBRARIES ${CRB_LIB_TOOLBOXMOR} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_toolboxmor_LIBRARIES PARENT_SCOPE)

feelpp_add_application(toolbox-opusheat SRCS toolboxmor_main.cpp LINK_LIBRARIES ${CRB_LIB_TOOLBOXMOR} Feelpp::feelpp_toolbox_heat_lib)
target_include_directories(feelpp_mor_toolbox-opusheat PUBLIC ${CMAKE_SOURCE_DIR}/toolboxes ${CMAKE_BINARY_DIR}/toolboxes)
if (EXISTS ${FEELPP_INSTALL_DIR}/include/feelpp/toolboxes)
  MESSAGE(STATUS "Loading toolboxes hearders from ${FEELPP_INSTALL_DIR}/include/feelpp/toolboxes")
  target_include_directories(feelpp_mor_toolbox-opusheat PUBLIC ${FEELPP_INSTALL_DIR}/include/feelpp/toolboxes)
endif()

feelpp_add_testcase(toolboxes PREFIX feelpp_crb CATEGORY CRB)

