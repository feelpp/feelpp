add_custom_target(feelpp_mor_apps)

feelpp_add_application( onlinerun SRCS crbonlinerun.cpp LINK_LIBRARIES  Feelpp::feelpp_mor INSTALL  MAN onlinerun EXEC ONLINERUN )
add_dependencies(feelpp_mor_apps feelpp_mor_onlinerun)
if ( TARGET Feelpp::feelpp_toolbox_heat_lib )

  mor_add_library( toolbox_heat
    SRCS toolboxmor_heat.cpp
    LINK_LIBRARIES  Feelpp::feelpp_toolbox_heat_lib_2dP1G1 Feelpp::feelpp_toolbox_heat_lib_3dP1G1 Feelpp::feelpp_mor
    EXEC CRB_LIB_TOOLBOX_HEAT PLUGIN )
  target_compile_definitions(${CRB_LIB_TOOLBOX_HEAT} PRIVATE FEELPP_MOR_PLUGIN_NAME=toolbox_heat FEELPP_MOR_PLUGIN_LIBNAME=${CRB_LIB_TOOLBOX_HEAT})



  feelpp_add_application(toolbox_heat_offline SRCS toolboxmor_heat_offline.cpp
    LINK_LIBRARIES
    ${CRB_LIB_TOOLBOX_HEAT}
    #feelpp_toolbox_heat_lib_2dP1G1
    #feelpp_toolbox_heat_lib_3dP1G1
    #Feelpp::feelpp_mor
    INSTALL EXEC MOR_TOOLBOX_HEAT)
  add_dependencies(feelpp_mor_apps ${MOR_TOOLBOX_HEAT} )

endif()
if (0)
if ( TARGET Feelpp::feelpp_toolbox_coefficientformpdes_lib)    
    feelpp_add_application(cfpdes SRCS toolboxmor_cfpdes.cpp LINK_LIBRARIES Feelpp::feelpp_toolbox_coefficientformpdes_lib Feelpp::feelpp_mor INSTALL EXEC MOR_TOOLBOX_CFPDES)
    add_dependencies(feelpp_mor_apps ${MOR_TOOLBOX_CFPDES} )
endif()
endif()
feelpp_add_testcase(opusheat PREFIX mor CATEGORY MOR)
