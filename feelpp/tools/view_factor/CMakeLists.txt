feelpp_add_testcase( cases )

feelpp_add_test( view_factor_cube NO_MPI_TEST  CFG cases/view_factor_cube.cfg GEO cases/view_factor_cube.geo )

feelpp_add_test( view_factor_circle NO_MPI_TEST  CFG cases/view_factor_cylinder.cfg GEO cases/view_factor_cylinder.geo)

feelpp_add_test( view_factor_subset NO_MPI_TEST  CFG cases/view_factor_cylinder.cfg GEO cases/view_factor_cylinder.geo)

