### TEMPLATE.txt.tpl; coding: utf-8 ---

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Natural convection problem
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

feelpp_add_application( room2d_crb
                        SRCS main.cpp ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp
                        CFG room2d.cfg GEO room2d.geo
                        DEFS CONVECTION_PARAMETER_SPACE_DIMENSION=4
                        EXCLUDE_FROM_ALL
                        )                            
feelpp_add_application( room3d_crb
                        SRCS main.cpp ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp
                        CFG room3d.cfg GEO room3d.geo
                        DEFS CONVECTION_PARAMETER_SPACE_DIMENSION=4 CONVECTION_DIM=3
                        EXCLUDE_FROM_ALL
                        )                            
