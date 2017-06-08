
add_library( feelpp-applications-crb-opusheat SHARED opusheat.cpp )
target_link_libraries( feelpp-applications-crb-opusheat ${FEELPP_LIBRARIES} )
install(FILES opusheat.hpp DESTINATION include/feelpp/crb/opusheat COMPONENT Devel)
install(TARGETS feelpp-applications-crb-opusheat DESTINATION lib/ COMPONENT Libs)
set(FEELPP_CRB_opusheat_LIBRARIES feelpp-applications-crb-opusheat PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_opusheat_LIBRARIES PARENT_SCOPE)

crb_add_model(opusheat OpusHeat HDRS opusheat.hpp
  LINK_LIBRARIES feelpp-applications-crb-opusheat )

feelpp_add_testcase(opusheat PREFIX feelpp_crb CATEGORY CRB)
