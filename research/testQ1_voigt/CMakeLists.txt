# cmake_minimum_required(VERSION 3.21.0)
# set(CMAKE_C_COMPILER clang)
# set(CMAKE_CXX_COMPILER clang++)

project(Avnir_feelpp)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# if(NOT CMAKE_BUILD_TYPE)
# set(CMAKE_BUILD_TYPE Release CACHE STRING "Choix du type de build" FORCE)
# endif()

# set(Feel++_DIR "/usr/share/feelpp/feel/cmake/modules")
# find_package(Feel++ REQUIRED)

# message(STATUS "Feel++ found: ${FEELPP_FOUND}")
# message(STATUS "Feel++_DIR: ${Feel++_DIR}")

# add_executable(voigt Q1_voigt.cpp)
# add_executable(Q1mat Q1_formMat.cpp)
# add_executable(Q2mat Q2_formMat.cpp)

# target_link_libraries(voigt PRIVATE Feelpp::feelpp)
# target_link_libraries(Q1mat PRIVATE Feelpp::feelpp)
# target_link_libraries(Q2mat PRIVATE Feelpp::feelpp)
add_subdirectory( test_expr_Q1_mat )
add_subdirectory( Q1_mat_with_voigt )