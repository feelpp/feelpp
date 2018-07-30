//+
h=0.1;
r=1.0;
L=2.0;
//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {0, -r, 0, h};
//+
Point(3) = {r, 0, 0, h};
//+
Point(4) = {L, 0, 0, h};
//+
Point(5) = {L, r, 0, h};
//+
Point(6) = {0, r, 0, h};
//+
Point(7) = {-r, 0, 0, h};
//+
Line(1) = {3, 4};
//+
Line(2) = {4, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 3};
//+
Circle(5) = {7, 1, 6};
//+
Circle(6) = {2, 1, 7};
//+
Circle(7) = {3, 1, 2};
//+
Line Loop(1) = {5, 4, 7, 6};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Line("Dirichlet") = {5, 6, 7, 4, 3, 1, 2};
//+
Physical Surface("Omega1") = {1};
//+
Physical Surface("Omega2") = {2};
