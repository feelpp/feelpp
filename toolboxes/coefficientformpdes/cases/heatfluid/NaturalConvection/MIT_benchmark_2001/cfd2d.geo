h=0.03;
Point(1) = {0., 0., 0., h};
Point(2) = {0, 8., 0., h};
Point(3) = {1, 8, 0., h};
Point(4) = {1., 0, 0., h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1,2,3, 4};
Plane Surface(1) = {1};

Physical Line("wall_cold") = {3};
Physical Line("insulated") = {2,4};
Physical Line("wall_hot") = {1};

Physical Surface("Omega") = {1};
