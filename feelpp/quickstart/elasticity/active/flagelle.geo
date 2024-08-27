h = 0.1;

Point(1) = {-0.2, 0, 0, h};
Point(2) = {0.2, 0, 0, h};
Point(3) = {0.2, 6.5, 0, h};
Point(4) = {-0.2, 6.5, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


Curve Loop(1) = {4, 1, 2, 3};
Surface(1) = {1};

Physical Curve("dirichlet", 5) = {1};
Physical Curve("Wall", 6) = {4, 3, 2};
Physical Surface("Caoutchouc", 7) = {1};