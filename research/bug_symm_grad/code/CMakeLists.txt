cmake_minimum_required(VERSION 3.21.0)
project( Bug_symm_grad )

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choix du type de build" FORCE)
endif()


add_executable(bug_symm_grad_tetra bug_symm_grad_tetra.cpp)
add_executable(bug_symm_grad_hexa bug_symm_grad_hexa.cpp)

target_link_libraries(bug_symm_grad_tetra PRIVATE Feelpp::feelpp)
target_link_libraries(bug_symm_grad_hexa PRIVATE Feelpp::feelpp)