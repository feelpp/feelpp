// Gmsh project created on Fri Jan  8 13:13:21 2021
//+
Point(1) = {0, 0, 0, 0.2};
//+
Point(2) = {1, 0, 0, 0.2};
//+
Point(3) = {1, 1, 0, 0.2};
//+
Point(4) = {0, 1, 0, 0.2};
//+
Point(5) = {0.2, 0.2, 0, 0.2};
//+
Point(6) = {0.2, 0.6, 0, 0.2};
//+
Point(7) = {0.6, 0.6, 0, 0.2};
//+
Point(8) = {0.6, 0.2, 0, 0.2};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 8};
//+
Line(6) = {8, 7};
//+
Line(7) = {7, 6};
//+
Line(8) = {6, 5};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {7, 8, 5, 6};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("FixedBoundaryMatTwo") = {6,7,8,5};

Physical Surface("Body") = {1};
