h = 0.25;
h = DefineNumber[ h, Name "Parameters/h" ];

L = 10;
H = 1;
B = 1;

Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, H, 0, h};
Point(4) = {0, H, 0, h};
Point(5) = {0, 0, B, h};
Point(6) = {L, 0, B, h};
Point(7) = {L, H, B, h};
Point(8) = {0, H, B, h};
Point(9) = {L, H/2, B/2, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Curve Loop(3) = {1, 10, -5, -9};
Plane Surface(3) = {3};
Curve Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};
Curve Loop(5) = {3, 12, -7, -11};
Plane Surface(5) = {5};
Curve Loop(6) = {4, 9, -8, -12};
Plane Surface(6) = {6};

Point{9} In Surface{4};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Physical Surface("Dirichlet") = {6};
Physical Point("tip") = {9};
Physical Volume("beam") = {1};
