add_custom_target(toolbox-mor)

crb_add_library( toolbox-mor SRCS toolboxmor.cpp EXEC CRB_LIB_TOOLBOXMOR )

set(FEELPP_CRB_toolboxmor_LIBRARIES ${CRB_LIB_TOOLBOXMOR} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_toolboxmor_LIBRARIES PARENT_SCOPE)

feelpp_add_application(toolboxmor-heat SRCS toolboxmor_heat.cpp LINK_LIBRARIES ${CRB_LIB_TOOLBOXMOR} Feelpp::feelpp_toolbox_heat_lib INSTALL EXEC MOR_TOOLBOX_HEAT)
feelpp_add_application(toolboxmor-cfpdes SRCS toolboxmor_cfpdes.cpp LINK_LIBRARIES ${CRB_LIB_TOOLBOXMOR} Feelpp::feelpp_toolbox_coefficientformpdes_lib INSTALL EXEC MOR_TOOLBOX_CFPDES)

feelpp_add_testcase(opusheat PREFIX feelpp_crb CATEGORY CRB)

add_dependencies(toolbox-mor ${MOR_TOOLBOX_HEAT} ${MOR_TOOLBOX_CFPDES})