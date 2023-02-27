h=0.05;
//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 5, 5, 0};
//
Rectangle(2) = {1, 1,0, 3, 3, 0};
//+
Curve Loop(3) = {3, 4, 1, 2};
//+
Curve Loop(4) = {7, 8, 5, 6};
//+
Plane Surface(3) = {3, 4};
//+
Physical Curve("OuterBdry") = {4, 3, 2, 1};
//+
Physical Curve("InnerTop") = {7};
//+
Physical Curve("LowerTop") = {5};
//+
Physical Curve("LeftTop") = {8};
//+
Physical Curve("RightTop") = {6};
//+
Physical Surface("Filling") = {3};

Characteristic Length { PointsOf{Line{5:8};}} = 0.5*h;//+

Characteristic Length { PointsOf{Line{1:4};}} = 4*h;//+
