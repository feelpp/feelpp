//+
h = 0.08;
SetFactory("OpenCASCADE");
Box(1) = {0, 0,0, 1.4, 1.4, 1.4};

Sphere(2) = {0.7, 0.7, 0.7, 0.08};



Delete{Volume{1};}

Surface Loop(4) = {1, 3, 5, 4, 2, 6};
Surface Loop(5) = {7};

Volume(4) = {4, 5};

Characteristic Length{ PointsOf{ Volume{2}; } } = 0.08;//+
Characteristic Length{ PointsOf{ Volume{4}; } } = 0.08;//+

Physical Surface("BoxWalls") = {3, 1, 4, 2,5,6};
Physical Surface("Ellipsoid") = {7};

Physical Volume("Fluid") = {4};
Physical Volume("EllipsoidVolume") = {2};