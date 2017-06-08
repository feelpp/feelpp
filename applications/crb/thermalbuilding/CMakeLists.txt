
add_library( feelpp-applications-crb-thermalbuilding SHARED thermalbuilding.cpp )
target_link_libraries( feelpp-applications-crb-thermalbuilding ${FEELPP_LIBRARIES} )
install(FILES thermalbuilding.hpp DESTINATION include/feelpp/crb/thermalbuilding COMPONENT Devel)
install(TARGETS feelpp-applications-crb-thermalbuilding DESTINATION lib/ COMPONENT Libs)
set(FEELPP_CRB_thermalbuilding_LIBRARIES feelpp-applications-crb-thermalbuilding PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_thermalbuilding_LIBRARIES PARENT_SCOPE)

crb_add_model(thermalbuilding ThermalBuilding HDRS thermalbuilding.hpp
  LINK_LIBRARIES feelpp-applications-crb-thermalbuilding )

feelpp_add_testcase(thermalbuilding PREFIX feelpp_crb CATEGORY CRB)
