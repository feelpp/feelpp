h = 0.4;
Point(1) = {0, 24, 0, h};
Point(2) = {0, 44, 0, h};
Point(3) = {0, 4, 0, h};

Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 3};
Curve Loop(1) = {1, 2};
Plane Surface(1) = {1};

Physical Curve("Wall", 4) = {1,2};
Physical Surface("Caoutchouc", 6) = {1};