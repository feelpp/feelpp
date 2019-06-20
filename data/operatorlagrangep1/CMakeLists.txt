
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

install( FILES ${OPLAGP1_MESHES} DESTINATION share/feelpp/data/operatorlagrangep1 )
