cmake_minimum_required(VERSION 3.21.0)

project(Linear_elasticity_with_Feelpp)

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choix du type de build" FORCE)
endif()

add_executable(elastic elastic.cpp)
target_link_libraries(elastic PRIVATE Feelpp::feelpp)