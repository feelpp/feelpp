h = 0.4;

Point(1) = {-25, 0, 0, h};
Point(2) = {-25, -1, 0, h};
Point(3) = {25, -1, 0, h};
Point(4) = {25, 0, 0, h};

Line(1) = {3, 2};
Line(2) = {2, 1};
Line(3) = {1, 4};
Line(4) = {4, 3};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

Physical Curve("boundary", 5) = {4, 1, 2, 3};
Physical Surface("wall", 6) = {1};
