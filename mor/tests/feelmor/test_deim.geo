h=0.1;

Point(1) = {0, 0, 0, h};
Point(2) = {0, 1, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {1, 0, 0, h};

Circle(1) = {4, 1, 2};
Line(2) = {2, 1};
Line(3) = {1, 4};
Line(4) = {4, 3};
Line(5) = {3, 2};

Line Loop(6) = {2, 3, 1};
Plane Surface(7) = {6};
Line Loop(8) = {5, -1, 4};
Plane Surface(9) = {8};

Physical Line("Gamma1") = {2};
Physical Line("Gamma2") = {5, 4, 3};
Physical Surface("Omega1") = {7};
Physical Surface("Omega2") = {9};
