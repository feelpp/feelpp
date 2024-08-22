h = 0.4;
Point(1) = {0, 24, 0, h};
Point(2) = {0, 44, 0, h};
Point(3) = {0, 4, 0, h};

Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 3};
Curve Loop(1) = {1, 2};
Plane Surface(1) = {1};

Point{1} In Surface{1};

Physical Point("CM")={1};
Physical Curve("Wall", 4) = {1,2};
Physical Surface("Caoutchouc", 6) = {1};

Point(100) = {-25, 0, 0, h};
Point(200) = {-25, -1, 0, h};
Point(300) = {25, -1, 0, h};
Point(400) = {25, 0, 0, h};

Line(100) = {300, 200};
Line(200) = {200, 100};
Line(300) = {100, 400};
Line(400) = {400, 300};

Curve Loop(100) = {400, 100, 200, 300};
Plane Surface(100) = {100};

Physical Curve("Obs1", 500) = {400, 100, 200, 300};
Physical Surface("obstacle", 600) = {100};