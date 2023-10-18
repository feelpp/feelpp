lc=0.06;
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {2, 0, 0, lc};
//+
Point(3) = {2, 6, 0, lc};
//+
Point(4) = {0, 6, 0, lc};
//+
Point(5) = {1, 4, 0, lc};
//+
Point(6) = {1.125, 4, 0, lc/3};
//+
Point(7) = {0.875, 4, 0, lc/3};
//+
Point(8) = {1, 4.125, 0, lc/3};
//+
Point(9) = {1, 3.875, 0, lc/3};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Circle(5) = {7, 5, 8};
//+
Circle(6) = {8, 5, 6};
//+
Circle(7) = {6, 5, 9};
//+
Circle(8) = {9, 5, 7};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Curve("Top") = {3};
//+
Physical Curve("Walls") = {2, 4};
//+
Physical Curve("CylinderSurface") = {5, 6, 7, 8};
//+
Physical Curve("Bottom") = {1};
//+
Physical Surface("Fluid") = {1};
//+
Physical Surface("CylinderVolume") = {2};
