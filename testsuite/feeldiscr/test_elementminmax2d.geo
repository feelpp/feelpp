SetFactory("OpenCASCADE");

Rectangle(1) = {-1, -1, 0, 2, 2, 0};
Disk(2) = {0, 0, 0, 0.5, 0.5};

Delete{Surface{1};}
Curve Loop(3) = {5};
Plane Surface(3) = {3};
Curve Loop(4) = {4, 1, 2, 3};
Curve Loop(5) = {5};
Plane Surface(4) = {4, 5};

Characteristic Length{ PointsOf{ Surface{2}; } } = 0.1;
Characteristic Length{ PointsOf{ Surface{4}; } } = 0.1;

Physical Curve("BoxWalls", 6) = {4, 3, 2, 1};
Physical Curve("CircleBoundary", 7) = {5};
Physical Surface("Circle", 8) = {2};
Physical Surface("Fluid", 9) = {4};
