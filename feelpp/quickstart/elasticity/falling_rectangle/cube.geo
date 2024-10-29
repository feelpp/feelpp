SetFactory("OpenCASCADE");
h = 1;
tol = 0.95;


Point(1) = {1, 1, 1, h};
Point(2) = {1, 1, -1, h};
Point(3) = {1, -1, 1, h};
Point(4) = {1, -1, -1, h};
Point(5) = {-1, -1, -1, h};
Point(6) = {-1, 1, -1, h};
Point(7) = {-1, 1, 1, h};
Point(8) = {-1, -1, 1, h};

Line(1) = {8, 5};
Line(2) = {8, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {8, 7};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {4, 2};
Line(9) = {1, 3};
Line(10) = {1, 2};
Line(11) = {2, 6};
Line(12) = {1, 7};

Curve Loop(1) = {12, -5, 2, -9};
Surface(1) = {1};
Curve Loop(3) = {10, -8, -3, -9};
Surface(2) = {3};
Curve Loop(5) = {4, -1, 2, 3};
Surface(3) = {5};
Curve Loop(7) = {5, -7, -6, -1};
Surface(4) = {7};
Curve Loop(9) = {11, -6, -4, 8};
Surface(5) = {9};
Curve Loop(11) = {10, 11, 7, -12};
Surface(6) = {11};

Surface Loop(1) = {1, 6, 4, 3, 2, 5};
Volume(1) = {1};

Point(9) = {0, -tol, tol, h};
Point(10) = {tol, 0, -tol, h};
Point(11) = {-tol, tol, 0, h};
Point(12) = {0, 0, 0, h};

Point{9} In Volume{1};
Point{10} In Volume{1};
Point{11} In Volume{1};
Point{12} In Volume{1};

Physical Point("a0", 1) = {12};
Physical Point("a1", 2) = {3};
Physical Point("a2", 3) = {2};
Physical Point("a3", 4) = {7};
Physical Surface("Wall", 13) = {4, 5, 2, 1};
Physical Volume("Caoutchouc", 14) = {1};
Physical Surface("Neumann1", 15) = {3};
Physical Surface("Neumann2", 16) = {6};