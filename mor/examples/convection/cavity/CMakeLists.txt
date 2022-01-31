### TEMPLATE.txt.tpl; coding: utf-8 ---


mor_add_library( cavity2d SRCS ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp   EXEC CRB_LIB_CAVITY2D )

install( FILES ../convection.hpp DESTINATION include/feelpp/crb/convection/cavity COMPONENT Devel)
set( FEELPP_CRB_CAVITY2D_LIBRARIES ${CRB_LIB_CAVITY2D} PARENT_SCOPE)
set( FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_CAVITY2D_LIBRARIES PARENT_SCOPE)

feelpp_add_application( cavity2dapp
                        SRCS main.cpp
                        LINK_LIBRARIES ${CRB_LIB_CAVITY2D}
                        CFG cavity.cfg GEO cavity2d.geo EXCLUDE_FROM_ALL )                            

# 3D Case
mor_add_library( cavity3d SRCS ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp   EXEC CRB_LIB_CAVITY3D )
target_compile_definitions( ${CRB_LIB_CAVITY3D} PRIVATE CONVECTION_DIM=3 )
set( FEELPP_CRB_CAVITY3D_LIBRARIES ${CRB_LIB_CAVITY3D} PARENT_SCOPE)
set( FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_CAVITY3D_LIBRARIES PARENT_SCOPE)

feelpp_add_application( cavity3dapp
                        SRCS main.cpp ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp
                        #LINK_LIBRARIES ${CRB_LIB_CAVITY3D}
                        DEFS CONVECTION_DIM=3
                        CFG cavity.cfg GEO cavity3d.geo )
