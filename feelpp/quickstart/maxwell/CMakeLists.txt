# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later

set_directory_properties(PROPERTIES LABELS "quickstart;maxwell")

add_custom_target(qs_maxwell)
add_dependencies(quickstart qs_maxwell)
add_dependencies(check qs_maxwell)

# Definite Maxwell problem with either Nedelec family and runtime order.
foreach(dim 2 3)
  feelpp_add_application( curlcurl_${dim}d SRCS qs_curlcurl.cpp INCLUDE_IN_ALL
    DEFS FEELPP_DIM=${dim} EXEC QS_CURLCURL INSTALL NO_TEST )
  add_dependencies(qs_maxwell ${QS_CURLCURL})
  foreach(family ned1 ned2)
    foreach(order 0 1 2 3 4)
      set(test_name feelpp_qs_curlcurl_${dim}d_${family}_k${order})
      add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${QS_CURLCURL}>
        --family=${family} --order=${order} --solution=polynomial
        --gmsh.hsize=0.5 --ksp-type=preonly --pc-type=lu
        --error-max=1e-8 --no-export --directory=testsuite/maxwell/${test_name})
      set_tests_properties(${test_name} PROPERTIES LABELS "quickstart;maxwell;hcurl" TIMEOUT 600)
      feelpp_append_solver_runtime_env(${test_name})
    endforeach()
    if(dim EQUAL 3)
      foreach(sphere_order 0 1 2)
        if(sphere_order EQUAL 0)
          set(sphere_error_max 8.0)
        elseif(sphere_order EQUAL 1)
          set(sphere_error_max 1.0)
        else()
          set(sphere_error_max 0.5)
        endif()
        foreach(sphere_geometry 1 2)
          set(test_name feelpp_qs_curlcurl_3d_${family}_sphere_k${sphere_order}_g${sphere_geometry})
          add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${QS_CURLCURL}>
            --config-file ${CMAKE_CURRENT_SOURCE_DIR}/cases/sphere/maxwell-sphere.cfg
            --family=${family} --order=${sphere_order} --geometry.order=${sphere_geometry}
            --error-max=${sphere_error_max} --no-export
            --directory=testsuite/maxwell/${test_name})
          set_tests_properties(${test_name} PROPERTIES LABELS "quickstart;maxwell;hcurl;sphere" TIMEOUT 600)
          feelpp_append_solver_runtime_env(${test_name})
        endforeach()
      endforeach()
    endif()
  endforeach()
endforeach()
