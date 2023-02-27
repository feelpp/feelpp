//+
h = 0.1;
height=2;
radius=1;
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 0, 0, height, radius, 2*Pi};
//Cylinder(1) = {-1,-1,-0.8, 0, 0, height, radius, 2*Pi};
//+
Physical Surface("TopDisk") = {2};
//+
Physical Surface("BottomDisk") = {3};
//+
Physical Surface("LateralSurface") = {1};
//+
Physical Volume("Cylinder") = {1};

Characteristic Length{ PointsOf{ Surface{1}; } } = h;