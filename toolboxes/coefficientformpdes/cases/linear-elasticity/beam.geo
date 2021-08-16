// horizontal beam with both extremity fixed
h=0.01;
Point(1) = {0, 0, 0, h};
Point(2) = {10, 0, 0, h};
Point(3) = {10, 0.3, 0, h};
Point(4) = {0, 0.3, 0, h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(7) = {6};
Physical Line("Dirichlet") = {2,4};
Physical Surface("Omega") = {7};