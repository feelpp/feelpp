
add_library( feelpp-applications-crb-heat3d SHARED heat3d.cpp )
target_link_libraries( feelpp-applications-crb-heat3d ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} )
INSTALL(TARGETS feelpp-applications-crb-heat3d DESTINATION lib/ COMPONENT Libs)

crb_add_model(heat3d Heat3d HDRS heat3d.hpp
  LINK_LIBRARIES feelpp-applications-crb-heat3d )

feelpp_add_testcase(heat3d PREFIX feelpp_crb CATEGORY CRB)
