set_directory_properties(PROPERTIES LABEL testcore )
IF ( TBB_FOUND )
  SET( OTESTS tbb )
ENDIF()
foreach(TEST traits singleton factory debug context simget ${OTESTS} env removecomments feelio)

  feelpp_add_test( ${TEST} )

endforeach()

feelpp_add_test( range NO_FEELPP_LIBRARY )

add_executable(test_boostmpi test_boostmpi.cpp)
target_link_libraries(test_boostmpi feelpp ${FEELPP_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
add_dependencies(testsuite test_boostmpi)
add_test(NAME test_boostmpi COMMAND mpirun -np 2 test_boostmpi --log_level=message  )


# test scatter/gather: identified a bug in PETSc up to 3.5.1
if(FEELPP_HAS_PETSC)
feelpp_add_test(gatherscatter_boostmpi_petsc_datatype SRCS test_gatherscatter.cpp DEFS USE_BOOST_MPI=1 USE_PETSC=1 USE_DATATYPE=1)
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

add_custom_target(feelpp_test_gatherscatter)
add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_mpi)
#add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_mpi_petsc)
add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_mpi_petsc_datatype)
add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_boostmpi)
#add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_boostmpi_petsc)
add_dependencies(feelpp_test_gatherscatter feelpp_test_gatherscatter_boostmpi_petsc_datatype)

#feelpp_add_application(mpi SRCS test_mpi.cpp)
feelpp_add_test( fusion )

if ( FEELPP_HAS_NT2 )
  feelpp_add_test(nt2)
endif()

if ( FEELPP_HAS_HDF5 )
  feelpp_add_test(hdf5)
endif()

if ( FEELPP_HAS_HARTS OR FEELPP_HAS_OPENMP )
  feelpp_add_test(multithread)
endif()
