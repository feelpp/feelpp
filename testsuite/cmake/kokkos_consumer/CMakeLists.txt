cmake_minimum_required(VERSION 3.21)
project(FeelppKokkosConsumer LANGUAGES C CXX)

find_package(Feel++ CONFIG REQUIRED)
enable_testing()

add_executable(feelpp_host_consumer host.cpp)
target_link_libraries(feelpp_host_consumer PRIVATE Feelpp::feelpp)
add_test(NAME feelpp_host_consumer COMMAND feelpp_host_consumer)

if(FEELPP_HAS_KOKKOS)
  if(NOT TARGET Feelpp::feelpp_kokkos)
    message(FATAL_ERROR "The installed Kokkos integration target is missing")
  endif()
  set_property(TARGET Kokkos::kokkos APPEND PROPERTY
    INTERFACE_COMPILE_DEFINITIONS FEELPP_KOKKOS_DEPENDENCY_SENTINEL=1)
  add_executable(feelpp_kokkos_consumer kokkos.cpp)
  target_link_libraries(feelpp_kokkos_consumer PRIVATE Feelpp::feelpp_kokkos)
  add_test(NAME feelpp_kokkos_consumer COMMAND feelpp_kokkos_consumer)
elseif(TARGET Feelpp::feelpp_kokkos)
  message(FATAL_ERROR "An installation without Kokkos exported a Kokkos integration target")
endif()
