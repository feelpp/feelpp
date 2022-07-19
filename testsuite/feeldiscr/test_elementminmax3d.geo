SetFactory("OpenCASCADE");
Box(1) = {-1, -1, -1, 2, 2, 2};
Sphere(2) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};

Delete{Volume{1};}
Surface Loop(4) = {1, 3, 5, 4, 2, 6};
Surface Loop(5) = {7};
Volume(4) = {4, 5};

Characteristic Length{ PointsOf{ Volume{2}; } } = 0.125;//+
Characteristic Length{ PointsOf{ Volume{4}; } } = 0.125;//+

Physical Surface("BoxWalls") = {3, 1, 4, 2,5,6};
Physical Surface("SphereBoundary") = {7};
Physical Volume("Fluid") = {4};
Physical Volume("Sphere") = {2};