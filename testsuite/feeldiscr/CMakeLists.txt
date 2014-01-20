set_directory_properties(PROPERTIES LABEL testdiscr )

feelpp_add_test( aitken NO_MPI_TEST )
feelpp_add_test( bdf )
feelpp_add_test( bdf2 CFG test_bdf2.cfg )

feelpp_add_test( disc NO_MPI_TEST )
feelpp_add_test( element )
feelpp_add_test( elementminmax )
feelpp_add_test( element_serialize NO_MPI_TEST )
feelpp_add_test( evaluator )

feelpp_add_test( functionallinearfree )
feelpp_add_test( fspace_context )
feelpp_add_test( geomap )

feelpp_add_test( lambda )
feelpp_add_test( lift )
feelpp_add_test( lm )
feelpp_add_test( lowerdim_entity )
feelpp_add_test( mesh_transfinite CFG test_mesh_transfinite.geocfg GEO test_mesh_transfinite.geo )
feelpp_add_test( mixed )
feelpp_add_test( normal3d )
feelpp_add_test( on_inside GEO test_on_inside.geo)

feelpp_add_test( operatorlinearfree )
feelpp_add_test( project )
feelpp_add_test( rbspace )
feelpp_add_test( rbspacevector )
feelpp_add_test( spaces )
feelpp_add_test( stencil )
feelpp_add_test( submesh )
feelpp_add_test( trace )
feelpp_add_test( vect_comp )

# feelpp_add_test( hypercubeinterpolation )
#feelpp_add_test( mortar )

feelpp_add_test( wire_basket NO_MPI_TEST )
feelpp_add_test( mortar NO_MPI_TEST)

