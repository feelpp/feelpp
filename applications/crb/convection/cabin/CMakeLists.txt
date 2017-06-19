### TEMPLATE.txt.tpl; coding: utf-8 ---

crb_add_library( cabin_lib SRCS ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp   EXEC CRB_LIB_CABIN )

install( FILES ../convection.hpp DESTINATION include/feelpp/crb/convection/cabin COMPONENT Devel)
set( FEELPP_CRB_CABIN_LIBRARIES ${CRB_LIB_CABIN} PARENT_SCOPE)
set( FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_CABIN_LIBRARIES PARENT_SCOPE)

feelpp_add_application( cabin
                        SRCS main.cpp
                        #LINK_LIBRARIES ${CRB_LIB_CABIN}
                        SRCS main.cpp ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp
                        CFG cabin.cfg GEO cabin.geo )                            
