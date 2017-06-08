
add_library( feelpp-applications-crb-heat3d SHARED heat3d.cpp )
target_link_libraries( feelpp-applications-crb-heat3d ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} )
install(FILES heat3d.hpp DESTINATION include/feelpp/crb/heat3d COMPONENT Devel)
install(TARGETS feelpp-applications-crb-heat3d DESTINATION lib/ COMPONENT Libs)
set(FEELPP_CRB_HEAT3D_LIBRARIES feelpp-applications-crb-heat3d PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_HEAT3D_LIBRARIES PARENT_SCOPE)

crb_add_model(heat3d Heat3d HDRS heat3d.hpp
  LINK_LIBRARIES feelpp-applications-crb-heat3d )

feelpp_add_testcase(heat3d PREFIX feelpp_crb CATEGORY CRB)
