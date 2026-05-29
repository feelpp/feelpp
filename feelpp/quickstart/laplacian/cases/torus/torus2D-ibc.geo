h = 0.05;
ri = 1.0;
re = 2.0;
delta = Pi/2.0;

Mesh.ElementOrder = 1;

Point(1) = {0, 0, 0, h};
Point(2) = {ri, 0, 0, h};
Point(3) = {re, 0, 0, h};
Point(4) = {re*Cos(delta), re*Sin(delta), 0, h};
Point(5) = {ri*Cos(delta), ri*Sin(delta), 0, h};

Line(1) = {2, 3};
Circle(2) = {3, 1, 4};
Line(3) = {4, 5};
Circle(4) = {5, 1, 2};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface("Omega") = {1};
Physical Curve("Dirichlet") = {1};
Physical Curve("Ibc") = {3};
Physical Curve("Neumann") = {2, 4};
