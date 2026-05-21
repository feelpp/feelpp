cmake_minimum_required(VERSION 3.21.0)
project( Elasticity_equation_testing_expression )

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choix du type de build" FORCE)
endif()

add_executable(Q1 Q1.cpp)
add_executable(Q2 Q2.cpp)

target_link_libraries(Q1 PRIVATE Feelpp::feelpp)
target_link_libraries(Q2 PRIVATE Feelpp::feelpp)