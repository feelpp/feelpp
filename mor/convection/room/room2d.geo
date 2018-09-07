cl1 = 1;
h=0.025;
l=0.25;
Point(1) = {0, 0, 0, h};
Point(3) = {1.5+l, 0, 0, h};
Point(31) = {1.5+l, 0.10, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {1.5, 1, 0, h};
Point(6) = {0, 0.9, 0, h};
Point(7) = {1.5, 0.10, 0, h};
Point(8) = {0.4, 0, 0, h};
Point(9) = {0.4, 0.25, 0, h};
Point(10) = {0.5, 0.25, 0, h};
Point(11) = {0.5, 0, 0, h};
Point(12) = {0.7, 0, 0, h};
Point(13) = {0.7, 0.5, 0, h};
Point(14) = {0.9, 0.5, 0, h};
Point(15) = {0.9, 0, 0, h};

Line(1) = {4, 6};
Line(2) = {6, 1};
Line(3) = {1, 8};
Line(4) = {8, 9};
Line(5) = {9, 10};
Line(6) = {10, 11};
Line(7) = {11, 12};
Line(8) = {12, 13};
Line(9) = {13, 14};
Line(10) = {14, 15};
Line(11) = {15, 3};
Line(12) = {3, 31};
Line(121) = {31, 7};
Line(13) = {7, 5};
Line(14) = {5, 4};

Line Loop(16) = {14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,121, 13};
Plane Surface(16) = {16};

Physical Line("inlet") = {1};
Physical Line("wall") = {2, 3, 7, 11,121, 13, 14};
Physical Line("cpu1") = {4, 5, 6};
Physical Line("cpu2") = {8, 9, 10};
Physical Line("outlet") = {12};

Physical Surface("omega") = {16};

