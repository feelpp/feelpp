
install(FILES opusheat.hpp DESTINATION include/feelpp/crb/opusheat COMPONENT Devel)
mor_add_library( opusheat SRCS opusheat.cpp  EXEC CRB_LIB_OPUSHEAT )
mor_add_library( opusheat-stationary SRCS opusheat-stationary.cpp  EXEC CRB_LIB_OPUSHEAT_STATIONARY )

set(FEELPP_CRB_opusheat_LIBRARIES ${CRB_LIB_OPUSHEAT} ${CRB_LIB_OPUSHEAT_STATIONARY} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_opusheat_LIBRARIES PARENT_SCOPE)

mor_add_model(opusheat OpusHeat CLASS OpusHeat<false> HDRS opusheat.hpp
  LINK_LIBRARIES ${CRB_LIB_OPUSHEAT} )

mor_add_model(opusheat-stationary OpusHeat CLASS OpusHeat<true> HDRS opusheat.hpp
  LINK_LIBRARIES ${CRB_LIB_OPUSHEAT_STATIONARY} )

feelpp_add_testcase(opusheat PREFIX feelpp_crb CATEGORY CRB)
