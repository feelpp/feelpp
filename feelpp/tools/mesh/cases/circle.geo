//+
SetFactory("OpenCASCADE");
Circle(1) = {1.8, 0.5, 0, 0.2, 0, 2*Pi};
//+
Physical Curve("BodyBoundary", 33) = {1};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("Body", 34) = {1};

Characteristic Length{ PointsOf{ Surface{1}; } } = 0.01;//+