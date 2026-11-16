# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later

set_directory_properties(PROPERTIES LABELS "testsuite;testsuite-feelpp;kokkos")

add_executable(feelpp_test_kokkos_runtime test_kokkos_runtime.cpp)
target_link_libraries(feelpp_test_kokkos_runtime PRIVATE Feelpp::feelpp_kokkos_runtime)
if(TARGET Feelpp::feelpp_kokkos)
  target_link_libraries(feelpp_test_kokkos_runtime PRIVATE Feelpp::feelpp_kokkos)
endif()
foreach(mode off lazy owned external petsc before-mpi after-mpi)
  string(REPLACE "-" "_" case_name "${mode}")
  add_test(NAME feelpp_test_kokkos_runtime_${mode}-np-1
    COMMAND feelpp_test_kokkos_runtime --run_test=runtime_${case_name})
  set_tests_properties(feelpp_test_kokkos_runtime_${mode}-np-1 PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OMP_PROC_BIND=false" TIMEOUT 60)
endforeach()
if(TARGET Feelpp::feelpp_kokkos)
  add_executable(feelpp_test_kokkos_smoke test_kokkos_smoke.cpp)
  target_link_libraries(feelpp_test_kokkos_smoke PRIVATE Feelpp::feelpp_kokkos)
  add_test(NAME feelpp_test_kokkos_smoke-np-1 COMMAND feelpp_test_kokkos_smoke)
  set_tests_properties(feelpp_test_kokkos_smoke-np-1 PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OMP_PROC_BIND=false" TIMEOUT 60)
endif()
add_dependencies(testsuite-feelpp feelpp_test_kokkos_runtime)
add_executable(feelpp_test_kokkos_environment test_kokkos_environment.cpp)
target_link_libraries(feelpp_test_kokkos_environment PRIVATE Feelpp::feelpp_kokkos_runtime)
set(environment_modes off lazy)
if(TARGET Feelpp::feelpp_kokkos)
  target_link_libraries(feelpp_test_kokkos_environment PRIVATE Feelpp::feelpp_kokkos)
  list(APPEND environment_modes owned external)
endif()
foreach(mode IN LISTS environment_modes)
  add_test(NAME feelpp_test_kokkos_environment_${mode}-np-1
    COMMAND feelpp_test_kokkos_environment --run_test=environment_${mode})
  set_tests_properties(feelpp_test_kokkos_environment_${mode}-np-1 PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OMP_PROC_BIND=false;OPENBLAS_NUM_THREADS=1" TIMEOUT 120)
endforeach()
add_dependencies(testsuite-feelpp feelpp_test_kokkos_environment)
if(TARGET feelpp_test_kokkos_smoke)
  add_dependencies(testsuite-feelpp feelpp_test_kokkos_smoke)
endif()
