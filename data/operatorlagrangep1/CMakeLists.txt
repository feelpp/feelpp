
set( OPLAGP1_MESHES
hypercube-2d-lagrange-order2.msh
simplex-2d-lagrange-order2.msh
simplex-2d-lagrange-order3-fekete.msh
simplex-2d-lagrange-order3-equispaced.msh
simplex-2d-lagrange-order4-fekete.msh
simplex-2d-lagrange-order4-equispaced.msh
simplex-2d-lagrange-order5-fekete.msh
simplex-2d-lagrange-order5-equispaced.msh
simplex-3d-lagrange-order2.msh
simplex-3d-lagrange-order3-warpblend.msh
simplex-3d-lagrange-order3-equispaced.msh
)

# Copy operator lagrangep1 meshes to build directory for testing
file(MAKE_DIRECTORY ${FEELPP_BUILD_DATADIR}/operatorlagrangep1)
foreach(mesh_file ${OPLAGP1_MESHES})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${mesh_file} 
                 ${FEELPP_BUILD_DATADIR}/operatorlagrangep1/${mesh_file} 
                 COPYONLY)
endforeach()

install( FILES ${OPLAGP1_MESHES} DESTINATION share/feelpp/data/operatorlagrangep1 )
