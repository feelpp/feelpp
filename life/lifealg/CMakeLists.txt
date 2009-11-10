add_library( lifealg
  dummy.cpp datamap.cpp 
  vector.cpp vectorublas.cpp 
  backend.cpp backendgmm.cpp backend_adaptive_reuse_pc.cpp backendpetsc.cpp 
  graphcsr.cpp matrixgmm.cpp solvereigen.cpp 
  solvernonlinear.cpp 
  )

target_link_libraries (lifealg lifecore)
