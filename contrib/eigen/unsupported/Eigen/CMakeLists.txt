set(Eigen_HEADERS AdolcForward BVH IterativeSolvers MatrixFunctions MoreVectorization AutoDiff AlignedVector3 Polynomials
                  CholmodSupport FFT NonLinearOptimization SparseExtra SuperLUSupport UmfPackSupport IterativeSolvers
                  NumericalDiff Skyline TaucsSupport MPRealSupport
   )

install(FILES
  ${Eigen_HEADERS}
  DESTINATION include/feel/unsupported/Eigen COMPONENT Devel
  )

add_subdirectory(src)
