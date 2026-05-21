cmake_minimum_required(VERSION 3.21.0)
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)

project(Linear_elasticity_with_Feelpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choix du type de build" FORCE)
endif()

set(Feel++_DIR "/usr/share/feelpp/feel/cmake/modules")
find_package(Feel++ REQUIRED)

message(STATUS "Feel++ found: ${FEELPP_FOUND}")
message(STATUS "Feel++_DIR: ${Feel++_DIR}")

add_executable(elastic elastic.cpp)

target_link_libraries(elastic PRIVATE Feelpp::feelpp)
