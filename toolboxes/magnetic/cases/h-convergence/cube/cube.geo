SetFactory("OpenCASCADE");
h=0.1;
Box(1) = {-1, -1, -1, 2, 2, 2};
Physical Surface("Boundary") = {3, 5, 2, 4, 1, 6};
Physical Volume("Omega") = {1};
MeshSize {4, 3, 7, 8, 1, 5, 6, 2} = h;

