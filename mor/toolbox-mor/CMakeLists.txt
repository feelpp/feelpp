add_custom_target(toolbox-mor)

add_library( feelpp_crb_toolbox-mor SHARED toolboxmor.cpp )
set_target_properties(feelpp_crb_toolbox-mor PROPERTIES VERSION 1 SOVERSION 1)
target_link_libraries( feelpp_crb_toolbox-mor PUBLIC Feelpp::feelpp   )
INSTALL(TARGETS feelpp_crb_toolbox-mor EXPORT feelpp-mor-export-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

set(FEELPP_CRB_toolboxmor_LIBRARIES ${CRB_LIB_TOOLBOXMOR} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_toolboxmor_LIBRARIES PARENT_SCOPE)

feelpp_add_application(toolboxmor-heat SRCS toolboxmor_heat.cpp LINK_LIBRARIES ${CRB_LIB_TOOLBOXMOR} Feelpp::feelpp_toolbox_heat_lib INSTALL EXEC MOR_TOOLBOX_HEAT)
feelpp_add_application(toolboxmor-cfpdes SRCS toolboxmor_cfpdes.cpp LINK_LIBRARIES ${CRB_LIB_TOOLBOXMOR} Feelpp::feelpp_toolbox_coefficientformpdes_lib INSTALL EXEC MOR_TOOLBOX_CFPDES)

feelpp_add_testcase(opusheat PREFIX feelpp_crb CATEGORY CRB)

add_dependencies(toolbox-mor ${MOR_TOOLBOX_HEAT} ${MOR_TOOLBOX_CFPDES})