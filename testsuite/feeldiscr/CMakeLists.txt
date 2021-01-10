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

# comment out as it breaks compilation and buildkite testuite pipeline 
# feelpp_add_test( space_0d SKIP_TEST )


feelpp_add_test( hdiv CFG test_hdiv.cfg NO_MPI_TEST)
feelpp_add_test( hdiv_oneelt CFG test_hdiv_oneelt.cfg NO_MPI_TEST )
# feelpp_add_test( hdiv3D CFG test_hdiv3D.cfg)
feelpp_add_test( hdiv3D_oneelt CFG test_hdiv3D_oneelt.cfg NO_MPI_TEST )
feelpp_add_test( hcurl_oneelt CFG test_hcurl_oneelt.cfg NO_MPI_TEST )
feelpp_add_test( hcurl_assembly CFG test_hcurl_assembly.cfg GEO two_elements_mesh.geo eight_elements_mesh.msh NO_MPI_TEST)
feelpp_add_test( hcurl CFG test_hcurl.cfg )
feelpp_add_test( hcurl_lag CFG test_hcurl.cfg DEFS HL_DIM=3 HL_ORDER=1)
feelpp_add_test( hcurl3D_oneelt CFG test_hcurl3D_oneelt.cfg NO_MPI_TEST )

feelpp_add_test( continuity_h1 )
feelpp_add_test( continuity_hdiv NO_MPI_TEST )
feelpp_add_test( continuity_hcurl NO_MPI_TEST )

feelpp_add_test( functionals NO_MPI_TEST )

feelpp_add_test( lambda )
feelpp_add_test( lift )
feelpp_add_test( lm )
feelpp_add_test( lowerdim_entity )
feelpp_add_test( meshstructured )
feelpp_add_test( mesh_transfinite CFG test_mesh_transfinite.cfg GEO test_mesh_transfinite.geo )
feelpp_add_test( mixed )
feelpp_add_test( normal3d )
feelpp_add_test( on_inside GEO test_on_inside.geo)

feelpp_add_test( operatorlinearfree )
feelpp_add_test( project NO_MPI_TEST)
feelpp_add_test( projector CLI "--pc-type=gasm --gmsh.hsize=0.2" ) # avoid destructor issue with LU-mumps in //
feelpp_add_test( spaces )

feelpp_add_test( stencil )

feelpp_add_test( makemesh )
feelpp_add_test( submesh )
feelpp_add_test( trace )
feelpp_add_test( element_component )
feelpp_add_test( element_component_3d )
feelpp_add_test( matrixfield )

if ( FEELPP_HAS_FFTW )
  feelpp_add_test( multiscaleimage)
endif()

feelpp_add_test( matching NO_MPI_TEST )

feelpp_add_test( wire_basket CLI "--gmsh.hsize=0.075" )
feelpp_add_test( mortar CFG test_mortar.cfg NO_MPI_TEST)
feelpp_add_test( mortar2 NO_MPI_TEST)

feelpp_add_test( extended )

feelpp_add_test( nlSolveComposite CFG test_nlSolveComposite.cfg)
feelpp_add_test( project_prodspace )

if ( FEELPP_HAS_SLEPC )
  feelpp_add_test( eigs )
endif()

feelpp_add_test( twolaplaciansdistributed )
feelpp_add_test( masterworker_laplacian )

feelpp_add_test( markerhdf5_2D SRCS test_markerhdf5 DEFS CHECKH5_DIM=2 GEO markerhdf5_2D.geo)
feelpp_add_test( markerhdf5_3D SRCS test_markerhdf5 DEFS CHECKH5_DIM=3 GEO markerhdf5_3D.geo)

feelpp_add_test( upmarker_faces )
feelpp_add_test( listofrange )

feelpp_add_test( productspaces CLI "--pc-type=gasm --sub-pc-factor-mat-solver-package-type=umfpack" )
feelpp_add_test( space_cst )
# feelpp_add_test( element0D )
feelpp_add_test( space_on_range )

feelpp_add_test( idelements)

feelpp_add_test( twospacesmap )

feelpp_add_test( on_dofs SRCS test_on_dofs GEO square_mmarked.geo CFG test_on_dofs.cfg )
#! not functional yet
# feelpp_add_test( element0D )

if ( FEELPP_HAS_VTK )
  feelpp_add_test( operatorlagrangep1 CLI "--pc-type=gamg --gmsh.hsize=0.25" )  
endif()

feelpp_add_test( sensors CFG sensors.cfg GEO sensorsdesc.json )
