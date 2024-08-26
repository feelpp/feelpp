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

Point(401) = {-35, 22, 0, h};
Point(402) = {35, -27, 0, h};
Point(403) = {-35, 19, 0, h};
Point(404) = {35, -30, 0, h};

Line(401) = {401, 402};
Line(402) = {401, 403};
Line(403) = {403, 404};
Line(404) = {404, 402};

Curve Loop(100) = {402, 403, 404, -401};
Plane Surface(100) = {100};

Physical Curve("Obs1", 500) = {402, 403, 404, -401};
Physical Surface("obstacle", 600) = {100};//+//+
