Point(1) = {0, 0, 0, 1.0};
Point(2) = {2, 0, 0, 1.0};
Point(3) = {2, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {4, 1, 2, 3};
Surface(1) = {1};

Physical Curve("Gamma_D", 5) = {4};
Physical Curve("Gamma_N", 6) = {2};
Physical Curve("Gamma_lat", 7) = {3, 1};
Physical Surface("Omega", 8) = {1};
