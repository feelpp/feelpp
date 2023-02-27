h = 0.01;
RDisk = 0.125;

// Construction of circle 
Point(1) = {0.5, 0.139, 0., h};
Point(2) = {0.5+RDisk, 0.139, 0., h};
Point(3) = {0.5-RDisk, 0.139, 0., h};

Circle(1) = {3,1,2};
Circle(2) = {2,1,3};

Line Loop(1) = {1,2};
Plane Surface(1) = {1};

Physical Curve("Circle") = {1,2};
Physical Surface("Cir") = {1};

// Rectangle vertices
Point(7) = {0.,0,0,h};
Point(8) = {1,0,0,h};
Point(9) = {1,0.6,0,h};
Point(10) = {0,0.6,0,h};

// Rectangle lines
Line(5) = {7, 8};
Line(6) = {8, 9};
Line(7) = {9, 10};
Line(8) = {10, 7};

Line Loop(3) = {5,6,7,8};
Plane Surface(3) = {3,1};

Characteristic Length{ PointsOf{ Surface{1}; } } = h;

Physical Curve("BoxWalls") = {8,7,6,5};
Physical Surface("Fluid") = {3};