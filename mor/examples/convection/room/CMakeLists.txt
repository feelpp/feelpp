### TEMPLATE.txt.tpl; coding: utf-8 ---

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Natural convection problem
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

mor_add_library( room2d SRCS ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp   EXEC CRB_LIB_ROOM2D )

install( FILES ../convection.hpp DESTINATION include/feelpp/crb/convection/room COMPONENT Devel)
set( FEELPP_CRB_ROOM2D_LIBRARIES ${CRB_LIB_ROOM2D} PARENT_SCOPE)
set( FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_ROOM2D_LIBRARIES PARENT_SCOPE)


feelpp_add_application( room2d_crb
                        SRCS main.cpp
                        CFG room2d.cfg GEO room2d.geo
                        DEFS CONVECTION_DIM=2
                        LINK_LIBRARIES ${CRB_LIB_ROOM2D}
                        )

add_compile_definitions(CONVECTION_DIM=3)
mor_add_library( room3d SRCS ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp   EXEC CRB_LIB_ROOM3D )
install( FILES ../convection.hpp DESTINATION include/feelpp/crb/convection/room COMPONENT Devel)
set( FEELPP_CRB_ROOM3D_LIBRARIES ${CRB_LIB_ROOM3D} PARENT_SCOPE )
set( FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_ROOM3D_LIBRARIES PARENT_SCOPE)

feelpp_add_application( room3d_crb
                        SRCS main.cpp
                        CFG room3d.cfg GEO room3d.geo
                        LINK_LIBRARIES ${CRB_LIB_ROOM3D}
                        DEFS CONVECTION_DIM=3
                        )
