//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-3, -3, 0, 6, 6, 0};

//+
Physical Surface("Fluid", 1) = {1};
Characteristic Length{ PointsOf{ Surface{1}; } } = 0.1;//+
Physical Curve("wall", 2) = {1,2,3,4};
