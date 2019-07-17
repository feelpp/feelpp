h = 0.1;

Point(1) = {0, 0, 0, h};
Point(2) = {0, 0, 5, h};

Line(1) = {1, 2};

Physical Point("inlet") = {1};
Physical Point("outlet") = {2};
Physical Line("centerline") = {1};
