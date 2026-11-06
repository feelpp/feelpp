h = 0.25;
h = DefineNumber[ h, Name "Parameters/h" ];

Point(1) = {0, 0, 0, h};
Point(2) = {10, 0, 0, h};
Point(3) = {10, 1, 0, h};
Point(4) = {0, 1, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Line("Dirichlet") = {4};
Physical Point("tip") = {3};
Physical Surface("beam") = {1};
