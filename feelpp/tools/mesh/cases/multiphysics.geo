h = 0.01;
SetFactory("OpenCASCADE");

Point(1) = {-1, 1, 0, h};
Point(2) = {-1, 0, 0, h};
Point(3) = {3, 1, 0, h};
Point(4) = {3, 0, 0, h};

Line(1) = {2, 4};
Line(2) = {4, 3};
Line(3) = {3, 1};
Line(4) = {1, 2};

Circle(5) = {-0.3, 0.5, 0, 0.2, 0, 2*Pi};
Circle(6) = {0.7, 0.5, 0, 0.2, 0, 2*Pi};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {6};
Curve Loop(3) = {5};
Plane Surface(1) = {1, 2, 3};
Curve Loop(4) = {5};
Plane Surface(2) = {4};
Curve Loop(5) = {6};
Plane Surface(3) = {5};

Physical Curve("BodyBoundary1", 7) = {5};
Physical Curve("BodyBoundary2", 8) = {6};
Physical Curve("BoxWalls", 9) = {1, 2, 3, 4};
Physical Surface("Fluid", 10) = {1};
Physical Surface("Body1", 11) = {2};
Physical Surface("Body2", 12) = {3};

Characteristic Length{ PointsOf{ Surface{2,3}; } } = 0.01;//+