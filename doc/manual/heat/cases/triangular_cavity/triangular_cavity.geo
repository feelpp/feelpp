h=0.025;

//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {4, 0, 0, h};
//+
Point(3) = {4, 1, 0, h};
//+
Point(4) = {0, 1, 0, h};
//+
Point(5) = {5, 1, 0, h};
//+
Point(6) = {5, 4, 0, h};
//+
Point(7) = {4, 4, 0, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {3, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 3};
//+
Line(9) = {4, 7};
//+
Point(9) = {-0.6, 1.8, 0, h};
//+
Point(10) = {3.4, 4.8, 0, h};
//+
Line(10) = {9, 4};
//+
Line(11) = {9, 10};
//+
Line(12) = {10, 7};
//+
Curve Loop(1) = {11, 12, -9, -10};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, 5, 6, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 4, 1, 2};
//+
Plane Surface(3) = {3};
//+
Physical Curve("RadiativeSurface1") = {3};
//+
Physical Curve("RadiativeSurface3") = {8};
//+
Physical Curve("RadiativeSurface2") = {9};
//+
Physical Surface("Rectangle2") = {1};
//+
Physical Surface("Rectangle3") = {2};
//+
Physical Surface("Rectangle1") = {3};
//+
Physical Curve("Fixed_T") = {1};
//+
Physical Curve("FixedQ2") = {11};
//+
Physical Curve("FixedQ3") = {6};