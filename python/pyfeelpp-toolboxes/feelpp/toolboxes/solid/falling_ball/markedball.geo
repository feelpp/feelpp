h = 0.4;
Point(1) = {0, 0, 0, h};
Point(2) = {20, 0, 0, h};
Point(3) = {-20, 0, 0, h};

Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 3};
Curve Loop(1) = {1, 2};
Plane Surface(1) = {1};



Physical Surface("Ball", 5) = {1};//+
Physical Curve("Gamma_C", 6) = {1};
//+
Physical Curve("Wall", 7) = {2};