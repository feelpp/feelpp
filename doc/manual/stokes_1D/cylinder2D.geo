h = 0.1;

Point(1) = {1, 0, 0, h};
Point(2) = {1, 0, 5, h};
Point(3) = {0, 0, 0, h};
Point(4) = {0, 0, 5, h};
Point(5) = {-1, 0, 0, h};
Point(6) = {-1, 0, 5, h};

Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {5, 6};
Line(4) = {1, 3, 5};
Line(5) = {2, 4, 6};

Line Loop(6) = {1, 5, -3, -4};
Plane Surface(6) = {6};

Physical Line("inlet") = {4};
Physical Line("wall") = {1, 3};
Physical Line("outlet") = {5};
Physical Surface("fluid") = {6};
Physical Line("centerline") = {2};
