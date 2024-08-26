h = 0.4;

Point(10) = {-29, -29, 0, h};
Point(20) = {-29, -24, 0, h};
Point(30) = {29, -24, 0, h};
Point(40) = {29, -29, 0, h};

Line(10) = {10, 40};
Line(20) = {40, 30};
Line(30) = {30, 20};
Line(40) = {20, 10};

Curve Loop(10) = {10, 20, 30, 40};
Plane Surface(10) = {10};

Physical Curve("Walls", 5) = {10, 20, 30, 40};
Physical Surface("Wall", 6) = {10};

