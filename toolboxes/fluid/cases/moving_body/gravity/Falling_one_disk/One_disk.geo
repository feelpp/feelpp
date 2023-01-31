h = 0.1;
lcCircle = 0.01;
lcDom = 0.1;
RDisk = 0.125;

// Construction of first disk

Point(1) = {1., 1.5, 0., lcCircle};
Point(2) = {1.+RDisk, 1.5, 0., lcCircle};
Point(3) = {1.-RDisk, 1.5, 0., lcCircle};

Circle(1) = {3,1,2};
Circle(2) = {2,1,3};

Line Loop(1) = {1,2};
Plane Surface(1) = {1};

Physical Curve("DiskFirst") = {1,2};
Physical Surface("DFirst") = {1};

// Rectangle vertices

Point(7) = {0,0,0,lcDom};
Point(8) = {2,0,0,lcDom};
Point(9) = {2,6,0,lcDom};
Point(10) = {0,6,0,lcDom};

// Rectangle lines

Line(5) = {7, 8};
Line(6) = {8, 9};
Line(7) = {9, 10};
Line(8) = {10, 7};

Line Loop(3) = {5,6,7,8};
Plane Surface(3) = {3,1};

Characteristic Length{ PointsOf{ Surface{1}; } } = lcCircle;

Physical Curve("BoxWalls") = {8,7,6,5};
Physical Surface("Fluid") = {3};