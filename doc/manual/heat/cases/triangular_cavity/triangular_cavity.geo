// Gmsh project created on Fri Dec 16 10:05:59 2022

h=0.025;

SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 4, 1, 0};
//+
Rectangle(2) = {3.99, 0.98, 0, 1, 3, 0};
//+
//+
Rectangle(3) = {0, 1, 0, 5, 1, 0};
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 36.87/180*Pi} {
  Surface{3}; 
}
//+
Translate {0.61, 0.19, 0} {
  Surface{3}; 
}
BooleanFragments{ Surface{3}; Delete; }{ Surface{2}; Delete; }
//+
BooleanFragments{ Surface{4}; Delete; }{ Surface{1}; Delete; }
//+
BooleanFragments{ Surface{6}; Delete; }{ Surface{2}; Delete; }

//+
Curve Loop(9) = {32, 12, 31, 30, 29, 33, -28};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {18, -17, -16, 15, 14, 13, -12};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {24, 28, 21, 27, 26, 25, 16, 17};
//+
Plane Surface(11) = {11};
//+
Physical Curve("RadiativeSurface1") = {24};
//+
Physical Curve("RadiativeSurface3") = {18};
//+
Physical Curve("RadiativeSurface2") = {32};
//+
Physical Surface("Rectangle2") = {9};
//+
Physical Surface("Rectangle3") = {10};
//+
Physical Surface("Rectangle1") = {11};
//+
Physical Curve("FixedQ3") = {29};
//+
Physical Curve("Fixed_T") = {26};
//+
Physical Curve("FixedQ2") = {14};

Delete{Surface{3:8};}

Characteristic Length{PointsOf{Line{:};}} = h;
