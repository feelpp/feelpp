h=0.1;
L=1;
H=4;
Point(1) = {0., 0., 0., h};
Point(2) = {L, 0., 0., h};
Point(3) = {L, H, 0., h};
Point(4) = {0., H, 0., h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};

Physical Line("outlet") = {3};
Physical Line("inlet") = {1};
Physical Line("left") = {4};
Physical Line("right") = {2};
Physical Surface("Omega") = {1};
