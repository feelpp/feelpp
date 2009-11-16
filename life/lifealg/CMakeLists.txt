if( PETSC_FOUND )
  set( LIFE_PETSC_SRCS vectorpetsc.cpp matrixpetsc.cpp solverlinearpetsc.cpp solvernonlinearpetsc.cpp)

  if ( HAVE_SLEPC )
    set( LIFE_SLEPC_SRCS solvereigenslepc.cpp )
  endif()
endif()

if ( TRILINOS_LIBRARIES )
  set( LIFE_TRILINOS_SRCS vectorepetra.cpp matrixepetra.cpp preconditionerml.cpp preconditionerifpack.cpp backendtrilinos.cpp solvernonlineartrilinos.cpp )
endif()

add_library( lifealg
  dummy.cpp datamap.cpp
  vector.cpp vectorublas.cpp
  backend.cpp backendgmm.cpp backend_adaptive_reuse_pc.cpp backendpetsc.cpp
  graphcsr.cpp matrixgmm.cpp solvereigen.cpp
  solvernonlinear.cpp
  ${LIFE_PETSC_SRCS}
  ${LIFE_SLEPC_SRCS}
  ${LIFE_TRILINOS_SRCS}
  )

target_link_libraries (lifealg lifecore)
