h = 0.01;

Point(1) = {-0.5, -0.5, 0, h};
Point(2) = {-0.5, 0.5, 0, h};
Point(3) = {0.5, 0.5, 0, h};
Point(4) = {0.5, -0.5, 0, h};
Point(5) = {0,0,0,h};

Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};

Curve Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};

Physical Curve("Wall", 8) = {1, 2, 3, 4};
Physical Surface("Caoutchouc", 9) = {1};
Physical Curve("Upper", 10) = {4};
Physical Curve("Lower", 11) = {2};
