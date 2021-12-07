add_custom_target(morapps)

set(FEELPP_CRB_toolboxmor_LIBRARIES feelpp_crb_toolbox-mor PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_toolboxmor_LIBRARIES PARENT_SCOPE)

add_subdirectory(onlinerun)

feelpp_add_application(toolboxmor-heat SRCS toolboxmor_heat.cpp LINK_LIBRARIES Feelpp::feelpp_mor Feelpp::feelpp_toolbox_heat_lib INSTALL EXEC MOR_TOOLBOX_HEAT)
feelpp_add_application(toolboxmor-cfpdes SRCS toolboxmor_cfpdes.cpp LINK_LIBRARIES Feelpp::feelpp_mor Feelpp::feelpp_toolbox_coefficientformpdes_lib INSTALL EXEC MOR_TOOLBOX_CFPDES)

feelpp_add_testcase(opusheat PREFIX feelpp_mor CATEGORY MOR)

add_dependencies(morapps ${MOR_TOOLBOX_HEAT} ${MOR_TOOLBOX_CFPDES})