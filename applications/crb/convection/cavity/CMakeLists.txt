### TEMPLATE.txt.tpl; coding: utf-8 ---

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Natural convection problem
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )



feelpp_add_application( cavity2d_crb
                        SRCS main.cpp ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp
                        CFG cavity.cfg GEO cavity2d.geo EXCLUDE_FROM_ALL
                        )                            

feelpp_add_application( cavity3d_crb
                        SRCS main.cpp ../residual.cpp ../jacobian.cpp ../other.cpp init.cpp specific.cpp
                        CFG cavity.cfg GEO cavity3d.geo
                        DEFS CONVECTION_DIM=3
                         EXCLUDE_FROM_ALL
                        )                            
