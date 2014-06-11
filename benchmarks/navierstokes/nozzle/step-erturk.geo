hmin=0.1;
hin = 1;
h = 1;
L=300*h;
Lin=20*h;
H=hin+h;
Point(1) =  { 0,       0, 0, hmin};
Point(2) =  { L,       0, 0, hmin};
Point(3) =  { L,       H, 0, hmin};
Point(4) =  { 0,       H, 0, hmin};
Point(5) =  { -Lin ,    H, 0, hmin};
Point(6) =  { -Lin ,    hin, 0, hmin};
Point(7) =  { 0 ,    hin, 0, hmin};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};

Line Loop(7) = {1, 2, 3, 4, 5, 6, 7};

Plane Surface(8) = {7};

Physical Line("inlet") = {5};
Physical Line("outlet") = {2};
Physical Line("wall") = {1,3, 4, 6,7};

Physical Surface(9) = {8};
