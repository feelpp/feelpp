set_directory_properties(PROPERTIES LABEL testmesh )
foreach(THETEST entity mesh regiontree mesh_codim1 kdtree P1mesh updatemarker partitioner_metis elementswithmarkedfaces meshmover convex meshfilters )

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

feelpp_add_test( meshmarker CFG test_meshmarker.cfg GEO test_meshmarker_square_4mat.geo )
feelpp_add_test( mesh_filters_operations SRCS test_mesh_filters_operations.cpp CFG test_mesh_filters_operations.cfg GEO test_mesh_filters_operations.geo )
feelpp_add_test( submesh1d SRCS test_submesh1d.cpp CFG test_submesh1d.cfg GEO test_submesh1d.geo )

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTS} )

if ( FEELPP_HAS_PARMMG OR FEELPP_HAS_MMG )
  feelpp_add_test( remesh CFG test_remesh.cfg  
                   GEO domains_2d.geo domains_3d.geo domains_2d_curved.geo domains_3d_curved.geo domain_2d.geo  domain_3d.geo 
                       test_remesh.json test_remesh_relations.json test_remesh_curved.json 
                   TIMEOUT 500)
endif()

feelpp_add_test(bvh2D SRCS test_bvh.cpp DEFS FEELPP_DIM=2 CFG test_bvh2d.cfg GEO domain_2d.geo NO_MPI_TEST )
feelpp_add_test(bvh3D SRCS test_bvh.cpp DEFS FEELPP_DIM=3 CFG test_bvh3d.cfg GEO cubic_cavity.geo NO_MPI_TEST )