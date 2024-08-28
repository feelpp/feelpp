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

Point(5) = {-0.2-1e-7, 10, 0, h};
Point(6) = {-0.2-1e-7, -3, 0, h};
Point(7) = {-0.5-1e-7, -3, 0, h};
Point(8) = {-0.5-1e-7, 10, 0, h};

Line(5) = {5, 8};
Line(6) = {8, 7};
Line(7) = {7, 6};
Line(8) = {6, 5};

Curve Loop(2) = {6, 7, 8, 5};
Surface(2) = {2};

Physical Surface("obstacle", 9) = {2};
Physical Curve("Obs1", 10) = {6, 5, 8, 7};