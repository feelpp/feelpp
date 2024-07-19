# feelpp_add_application( grepldeim
#                         SRCS main.cpp
#                         CFG grepldeim.cfg
#                         LINK_LIBRARIES Feelpp::feelpp_mor
#                         )

mor_add_library( grepldeim SRCS main.cpp EXEC CRB_LIB_GREPLDEIM PLUGIN )
target_compile_definitions(${CRB_LIB_GREPLDEIM} PRIVATE FEELPP_MOR_PLUGIN_NAME=grepldeim FEELPP_MOR_PLUGIN_LIBNAME=${CRB_LIB_GREPLDEIM})
install(FILES grepldeim.hpp DESTINATION include/feelpp/crb/grepldeim COMPONENT Devel)
set(FEELPP_CRB_grepldeim_LIBRARIES ${CRB_LIB_GREPLDEIM} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_grepldeim_LIBRARIES PARENT_SCOPE)

mor_add_model(grepldeim GreplDEIM CLASS GreplDEIM<2,2> HDRS grepldeim.hpp
  LINK_LIBRARIES ${CRB_LIB_GREPLDEIM} CFG grepldeim.cfg)
