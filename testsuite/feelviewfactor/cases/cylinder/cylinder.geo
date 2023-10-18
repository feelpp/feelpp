//+
h = 0.1;
height=2;
radius=1;
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, -0.2, 0, 0, height+0.4, radius+0.2,2*Pi};
Cylinder(2) = {0, 0, 0, 0, 0, height, radius, 2*Pi};
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

//+
Physical Surface("TopDisk") = {5};
//+
Physical Surface("BottomDisk") = {6};
//+
Physical Surface("LateralSurface") = {4};
//+
Physical Volume("Cylinder") = {1};

Characteristic Length{ PointsOf{ Volume{:}; } } = h;//+
