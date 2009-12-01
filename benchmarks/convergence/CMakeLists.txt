
INCLUDE_DIRECTORIES (
  ${LIFE_SOURCE_DIR}/benchmarks/convergence
  )
ENABLE_TESTING()

foreach( LAPLACIAN laplacian1dP1 laplacian1dP5
    laplacian1d2dP1 laplacian1d2dP5
    laplacian2dP1 laplacian2dP5
    laplacian2d3dP1
    laplacian2dP1h laplacian2dP5h
    laplacian3dP1 laplacian3dP2 laplacian3dP3 laplacian3dP5
    laplacian3dP1h
    laplacian2dP1v laplacian2dP4v )

    life_add_bench( ${LAPLACIAN} )

  add_test( ${LAPLACIAN} 
      ${LIFE_SOURCE_DIR}/benchmarks/python/batch_validation.py
      --program=${targetname}
      --path="${LIFE_BUILD_DIR}/benchmarks/convergence/"
      --cmdargs="--backend=petsc" )
endforeach()

foreach( bench stokes nonlinearpow )
   life_add_bench( ${bench} )
   add_test( ${bench}
      ${LIFE_SOURCE_DIR}/benchmarks/python/batch_validation.py
      --program=${targetname}
      --path="${LIFE_BUILD_DIR}/benchmarks/convergence/"
      --cmdargs="--backend=petsc" )
endforeach()

foreach( ALE ale1 ale2 ale3 ale4 ale5 )
  life_add_bench( ${ALE} )
  
  add_test( ${ALE}
      ${LIFE_SOURCE_DIR}/benchmarks/python/batch_validation.py
      --program=${targetname}
      --path="${LIFE_BUILD_DIR}/benchmarks/convergence/"
      --cmdargs="--backend=petsc" )
endforeach()
