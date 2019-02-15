
set( OPLAGP1_MESHES
hypercube-2d-lagrange-order2.msh
simplex-2d-lagrange-order2.msh
simplex-3d-lagrange-order2.msh
)

install( FILES ${OPLAGP1_MESHES} DESTINATION share/feelpp/data/operatorlagrangep1 )
