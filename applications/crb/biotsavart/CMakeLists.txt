include_directories(${CMAKE_SOURCE_DIR}/applications/crb/thermoelectric)

feelpp_add_application( biotsavartCRB
                        SRCS biotsavart_main.cpp biotsavart.hpp )
feelpp_add_application( biotsavartCRBOpt
                        SRCS biotsavart_opt.cpp biotsavart.hpp )

feelpp_add_testcase( te-crb-quarterturn )
