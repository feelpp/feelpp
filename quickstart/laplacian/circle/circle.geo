// Gmsh project created on Fri Mar 18 16:22:09 2016
h=0.1;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {-1, 0, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {0, -1, 0, h};
Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(5) = {5, 1, 2};
Line Loop(6) = {1, 2, 3, 5};
Plane Surface(7) = {6};
Physical Line("Dirichlet") = {1, 2};
Physical Line("Neumann") = {3};
Physical Line("Robin") = {5};
Physical Surface(8) = {7};
