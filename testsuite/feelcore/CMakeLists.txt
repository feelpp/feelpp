set_directory_properties(PROPERTIES LABEL testcore )
IF ( TBB_FOUND )
  SET( OTESTS tbb )
ENDIF()
foreach(TEST json git traits singleton enumerate factory debug context simget ${OTESTS} env removecomments feelio boostmpi enums mongo hana)

  feelpp_add_test( ${TEST} )

endforeach()
if ( SPECX_FOUND )
  feelpp_add_test(taskflow LINK_LIBRARIES Feelpp::feelpp_taskflow)
endif()
feelpp_add_test( range )

feelpp_get_compile_definition(Feelpp::feelpp FEELPP_HAS_PETSC)
# test scatter/gather: identified a bug in PETSc up to 3.5.1
if(FEELPP_HAS_PETSC)
feelpp_add_test(gatherscatter_boostmpi_petsc_datatype SRCS test_gatherscatter.cpp DEFS USE_BOOST_MPI=1 USE_PETSC=1 USE_DATATYPE=1  )
endif()
# the next one fails for now
#feelpp_add_test(gatherscatter_boostmpi_petsc SRCS test_gatherscatter.cpp DEFS USE_BOOST_MPI=1 USE_PETSC=1)
feelpp_add_test(gatherscatter_boostmpi  SRCS test_gatherscatter.cpp DEFS USE_BOOST_MPI=1)
if(FEELPP_HAS_PETSC)
feelpp_add_test(gatherscatter_mpi_petsc_datatype  SRCS test_gatherscatter.cpp DEFS USE_PETSC=1 USE_DATATYPE=1)
endif()
# the next one fails
#feelpp_add_test(gatherscatter_mpi_petsc  SRCS test_gatherscatter.cpp DEFS USE_PETSC=1)
feelpp_add_test(gatherscatter_mpi  SRCS test_gatherscatter.cpp )
feelpp_add_test(comm  SRCS test_comm.cpp )
feelpp_add_test(worldcomm  SRCS test_worldcomm.cpp GEO test_twodomains.geo )

#add_custom_target(feelpp_test_gatherscatter)
#add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_mpi)
#add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_mpi_petsc)
#add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_boostmpi)
#if(FEELPP_HAS_PETSC)
#add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_mpi_petsc_datatype)
#add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_boostmpi_petsc_datatype)
#endif()

#add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_boostmpi_petsc)


#feelpp_add_application(mpi SRCS test_mpi.cpp)
feelpp_add_test( fusion )

if ( FEELPP_HAS_NT2 )
  feelpp_add_test(nt2)
endif()

feelpp_get_compile_definition(Feelpp::feelpp FEELPP_HAS_HDF5)
if ( FEELPP_HAS_HDF5 )
  message(STATUS "compiling feelpp_test_hdf5...")
  feelpp_add_test(hdf5)
endif()

if ( FEELPP_HAS_HARTS OR FEELPP_HAS_OPENMP )
  feelpp_add_test(multithread)
endif()

feelpp_add_test(serialization)
feelpp_add_test(ptree_editions CFG test_ptree_editions.cfg GEO test_ptree_editions.json )

feelpp_add_test( observer )
feelpp_add_test( observer_mesh )
feelpp_add_test( eigenrand )
