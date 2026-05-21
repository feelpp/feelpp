cmake_minimum_required(VERSION 3.21.0)
project( Test_Q1_voigt )

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choix du type de build" FORCE)
endif()

add_executable(Q1voigt Q1_voigt.cpp)
add_executable(Q1mandel Q1_mandel.cpp)

target_link_libraries(Q1voigt PRIVATE Feelpp::feelpp)
target_link_libraries(Q1mandel PRIVATE Feelpp::feelpp)