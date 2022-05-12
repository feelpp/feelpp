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
//+
Extrude {0, 0, 1} {
  Curve{1}; Curve{2}; Curve{3}; Curve{4}; 
}
//+
Translate {0, 0, 0.2} {
  Curve{7}; Curve{6}; Curve{5}; Curve{8}; 
}
//+
Extrude {0, 0, 0.4} {
  Curve{7}; Curve{8}; Curve{5}; Curve{6}; 
}
//+
Curve Loop(1) = {13, -15, -2, 11};
//+
Curve Loop(2) = {17, 21, 9, 13};
//+
Plane Surface(41) = {2};
//+
Curve Loop(3) = {3, 4, 1, 2};
//+
Plane Surface(42) = {3};
//+
Curve Loop(4) = {25, 29, 33, 37};
//+
Plane Surface(43) = {4};
//+
Curve Loop(5) = {6, 7, 8, 5};
//+
Plane Surface(44) = {5};
//+
Surface Loop(1) = {20, 42, 24, 12, 16, 41};
//+
Surface Loop(2) = {40, 44, 28, 32, 36, 43};
//+
Volume(1) = {1, 2};
//+
Physical Surface("OtherFixedBoundary") = {43, 28, 32, 40, 44, 36};

Physical Volume("Body") = {1};