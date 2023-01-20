//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
//+
Physical Curve("BodyBoundary", 3) = {1};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("Body", 4) = {1};

Characteristic Length{ PointsOf{ Surface{1}; } } = 0.04;//+