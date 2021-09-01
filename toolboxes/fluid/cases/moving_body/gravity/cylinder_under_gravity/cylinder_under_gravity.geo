Point(1) = {0, 0, 0, 0.05};
//
Point(2) = {2, 0, 0, 0.05};
//+
Point(3) = {2, 6, 0, 0.05};
//+
Point(4) = {0, 6, 0, 0.05};
//+

//+
Point(5) = {1, 4, 0, 0.1};
//+
Point(6) = {1.125, 4, 0, 0.01};

//+
Point(8) = {0.875, 4, 0, 0.01};
//+
Point(9) = {1, 4.125, 0, 0.01};
//+
Point(10) = {1, 3.875, 0, 0.01};
//+
Circle(1) = {9, 5, 6};
//+
Circle(2) = {6, 5, 10};
//+
Circle(3) = {10, 5, 8};
//+
Circle(4) = {8, 5, 9};
//+
Line(5) = {4, 3};
//+
Line(6) = {3, 2};
//+
Line(7) = {2, 1};
//+
Line(8) = {1, 4};
//+
Curve Loop(1) = {5, 6, 7, 8};
//+
Curve Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Surface("Fluid") = {1};
//+
Physical Surface("Swimmer") = {2}; // Volume of the cylinder
//+
Physical Curve("Walls") = {8, 6};
//+
Physical Curve("Top") = {5};
//+
Physical Curve("Bottom") = {7};
Physical Curve("Head") = {4, 1, 2, 3}; // Boundary of the cylinder
