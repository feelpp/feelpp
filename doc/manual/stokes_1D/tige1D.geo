h = 0.1;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {2, 0, 0, h};
Point(4) = {5, 0, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

Physical Point("inlet") = {1};
Physical Point("outlet") = {4};
Physical Line("centerline") = {4};