set_directory_properties(PROPERTIES LABEL testdiscr )

feelpp_add_test( aitken NO_MPI_TEST )

feelpp_add_test( disc NO_MPI_TEST )
feelpp_add_test( element )
feelpp_add_test( elementminmax )
feelpp_add_test( element_serialize NO_MPI_TEST )
feelpp_add_test( evaluator )

feelpp_add_test( functionallinearfree )
feelpp_add_test( fspace_context )
feelpp_add_test( geomap )

feelpp_add_test( spacefunc )


feelpp_add_test( hdiv CFG test_hdiv.cfg NO_MPI_TEST)
feelpp_add_test( hdiv_oneelt CFG test_hdiv_oneelt.cfg NO_MPI_TEST )
feelpp_add_test( hdiv3D CFG test_hdiv3D.cfg)
feelpp_add_test( hdiv3D_oneelt CFG test_hdiv3D_oneelt.cfg NO_MPI_TEST )
feelpp_add_test( hcurl_oneelt CFG test_hcurl_oneelt.cfg NO_MPI_TEST )
feelpp_add_test( hcurl_assembly CFG test_hcurl_assembly.cfg GEO two_elements_mesh.geo eight_elements_mesh.msh NO_MPI_TEST)
feelpp_add_test( hcurl CFG test_hcurl.cfg )
feelpp_add_test( hcurl_lag CFG test_hcurl.cfg DEFS HL_DIM=3 HL_ORDER=1)
feelpp_add_test( hcurl3D_oneelt CFG test_hcurl3D_oneelt.cfg NO_MPI_TEST )

feelpp_add_test( continuity_h1 NO_MPI_TEST )
feelpp_add_test( continuity_hdiv NO_MPI_TEST )
feelpp_add_test( continuity_hcurl NO_MPI_TEST )

feelpp_add_test( functionals NO_MPI_TEST )

feelpp_add_test( lambda )
feelpp_add_test( lift )
feelpp_add_test( lm )
feelpp_add_test( lowerdim_entity )
feelpp_add_test( mesh_transfinite CFG test_mesh_transfinite.cfg GEO test_mesh_transfinite.geo )
feelpp_add_test( mixed )
feelpp_add_test( normal3d )
feelpp_add_test( on_inside GEO test_on_inside.geo)

feelpp_add_test( operatorlinearfree )
feelpp_add_test( project NO_MPI_TEST)
feelpp_add_test( spaces )
feelpp_add_test( stencil )
feelpp_add_test( submesh )
feelpp_add_test( trace )
feelpp_add_test( vect_comp )
feelpp_add_test( matching NO_MPI_TEST )

feelpp_add_test( wire_basket NO_MPI_TEST )
feelpp_add_test( mortar CFG test_mortar.cfg NO_MPI_TEST)

feelpp_add_test( exporter   CFG exporter.cfg NO_MPI_TEST)
feelpp_add_test( extended )

feelpp_add_test( nlSolveComposite CFG test_nlSolveComposite.cfg)
feelpp_add_test( project_prodspace )

if ( FEELPP_HAS_SLEPC )
  feelpp_add_test( eigs )
endif()
