hmin=0.00017;

Point(1) =  { 0,       0.002, 0, hmin};
Point(2) =  { 0,      -0.002, 0, hmin};
Point(3) =  { 0.008,  -0.002, 0, hmin};
Point(4) =  { 0.008,  -0.006, 0, hmin};
Point(5) =  { 0.12 ,  -0.006, 0, hmin};
Point(6) =  { 0.12 ,   0.002, 0, hmin};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line Loop(7) = {1, 2, 3, 4, 5, 6};

Plane Surface(8) = {7};

Physical Line("inlet") = {1};
Physical Line("outlet") = {5};
Physical Line("wall") = {2, 3, 4, 6};

Physical Surface(9) = {8};
