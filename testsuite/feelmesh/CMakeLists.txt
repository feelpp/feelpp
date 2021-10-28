set_directory_properties(PROPERTIES LABEL testmesh )
foreach(THETEST entity mesh regiontree mesh_codim1 kdtree P1mesh meshmarker updatemarker partitioner_metis elementswithmarkedfaces meshmover convex )

  if(THETEST MATCHES partitioner_metis)
    feelpp_get_compile_definition(Feelpp::feelpp_contrib FEELPP_HAS_METIS)
    if(FEELPP_HAS_METIS)
      message(STATUS "compiling feelpp_test_partitioner_metis...")
      feelpp_add_test( ${THETEST} )
    endif()
  else()
    feelpp_add_test( ${THETEST} )
  endif()

endforeach()

feelpp_add_test( mesh_filters_operations SRCS test_mesh_filters_operations.cpp CFG test_mesh_filters_operations.cfg GEO test_mesh_filters_operations.geo )
feelpp_add_test( submesh1d SRCS test_submesh1d.cpp CFG test_submesh1d.cfg GEO test_submesh1d.geo )

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )

if ( FEELPP_HAS_PARMMG OR FEELPP_HAS_MMG )
  feelpp_add_test( remesh CLI "--gmsh.hsize=0.2 --functions.s=0.1 --backend.rebuild=true --pc-type gamg --ksp-monitor=1 --pc-view=1" NO_MPI_TEST TIMEOUT 400)
endif()
