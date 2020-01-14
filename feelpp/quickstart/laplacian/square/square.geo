// Gmsh project created on Fri Mar 18 16:22:09 2016
h=0.1;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(4) = {1, 2, 3, 4};
Plane Surface(5) = {4};
Physical Line("Dirichlet") = {1,2,4};
Physical Line("Neumann") = {3};
Physical Surface("Omega") = {5};
