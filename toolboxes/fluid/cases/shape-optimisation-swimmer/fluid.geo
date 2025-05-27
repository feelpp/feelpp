//+
h = 0.08;
SetFactory("OpenCASCADE");
Box(1) = {0, 0,0, 1, 1, 1};

Sphere(2) = {0.5, 0.5, 0.5, 0.08};


Delete{Volume{1};}

Surface Loop(4) = {1, 3, 5, 4, 2, 6};
Surface Loop(5) = {7};

Volume(4) = {4, 5};

Delete{Volume{2};}

Characteristic Length{ PointsOf{ Volume{4}; } } = 0.08;//+

Physical Surface("BoxWalls") = {3, 1, 4, 2,5,6};
Physical Surface("Ellipsoid") = {7};

Physical Volume("Fluid") = {4};