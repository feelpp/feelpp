add_custom_target(feelpp_mor_apps)

feelpp_add_application( onlinerun SRCS crbonlinerun.cpp  INSTALL  MAN onlinerun EXEC ONLINERUN )
add_dependencies(feelpp_mor_apps feelpp_mor_onlinerun)

if ( TARGET Feelpp::feelpp_toolbox_heat_lib )
    feelpp_add_application(heat SRCS toolboxmor_heat.cpp LINK_LIBRARIES  Feelpp::feelpp_toolbox_heat_lib Feelpp::feelpp_mor INSTALL EXEC MOR_TOOLBOX_HEAT)
    add_dependencies(feelpp_mor_apps ${MOR_TOOLBOX_HEAT} )
endif()
if ( TARGET Feelpp::feelpp_toolbox_coefficientformpdes_lib)    
    feelpp_add_application(cfpdes SRCS toolboxmor_cfpdes.cpp LINK_LIBRARIES Feelpp::feelpp_toolbox_coefficientformpdes_lib Feelpp::feelpp_mor INSTALL EXEC MOR_TOOLBOX_CFPDES)
    add_dependencies(feelpp_mor_apps ${MOR_TOOLBOX_CFPDES} )
endif()

feelpp_add_testcase(opusheat PREFIX mor CATEGORY MOR)
