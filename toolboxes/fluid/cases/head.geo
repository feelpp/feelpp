h = 0.05;

l = 0.5;
ha = 1.5;

hbox = 0.2;
habox = 6;

Point(1) = {0 - l/2, 1 - ha/2, 0, h};
Point(2) = {0 - l/2, 1 + ha/2, 0, h};
Point(3) = {0 + l/2, 1 + ha/2, 0, h};
Point(4) = {0 + l/2, 1 - ha/2, 0, h};

Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};

Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};

Point(11) = {-1, 1 - habox/2, 0, hbox};
Point(12) = {8., 1 - habox/2, 0, hbox};
Point(13) = {-1, 1 + habox/2, 0, hbox};
Point(14) = {8., 1 + habox/2, 0, hbox};

Line(5) = {13, 11};
Line(6) = {11, 12};
Line(7) = {12, 14};
Line(8) = {14, 13};

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {1, 2};

Physical Curve("boundary_head", 8) = {4, 1, 2, 3};
Physical Surface("head", 11) = {1};
Physical Curve("BoxWalls", 12) = {5, 8, 7, 6};
Physical Surface("Fluid", 13) = {2};
