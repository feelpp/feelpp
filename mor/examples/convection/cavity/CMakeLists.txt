### TEMPLATE.txt.tpl; coding: utf-8 ---
mor_add_library( cavity2d SRCS ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp DEFS CONVECTION_DIM=2 EXEC CRB_LIB_CAVITY2D )
install( FILES ../convection.hpp DESTINATION include/feelpp/crb/convection/cavity COMPONENT Devel)
set( FEELPP_CRB_CAVITY2D_LIBRARIES ${CRB_LIB_CAVITY2D} PARENT_SCOPE)
set( FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_CAVITY2D_LIBRARIES PARENT_SCOPE)

feelpp_add_application( cavity2dapp
                        SRCS main.cpp
                        DEFS CONVECTION_DIM=2
                        CFG cavity.cfg GEO cavity2d.geo
                        LINK_LIBRARIES ${CRB_LIB_CAVITY2D})

# 3D Case
mor_add_library( cavity3d SRCS ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp DEFS CONVECTION_DIM=3 EXEC CRB_LIB_CAVITY3D )
set( FEELPP_CRB_CAVITY3D_LIBRARIES ${CRB_LIB_CAVITY3D} PARENT_SCOPE)
set( FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_CAVITY3D_LIBRARIES PARENT_SCOPE)

feelpp_add_application( cavity3dapp
                        SRCS main.cpp
                        DEFS CONVECTION_DIM=3
                        CFG cavity.cfg GEO cavity3d.geo
                        LINK_LIBRARIES ${CRB_LIB_CAVITY3D})
