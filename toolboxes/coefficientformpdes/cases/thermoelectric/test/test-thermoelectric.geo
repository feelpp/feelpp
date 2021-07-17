h=0.1;
Point(1) = {2.0000000000000001e-01,0.0000000000000000e+00,0.0000000000000000e+00,h};
Point(2) = {0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,h};
Point(3) = {-2.0000000000000001e-01,0.0000000000000000e+00,0.0000000000000000e+00,h};
Point(4) = {-1,-1,0,h};
Point(5) = {1,-1,0,h};
Point(6) = {1,1,0,h};
Point(7) = {-1,1,0,h};
thickness=0.4;
Point(8) = {-1,1+thickness,0,h};
Point(9) = {1,1+thickness,0,h};
Point(10) = {-1,-1-thickness,0,h};
Point(11) = {1,-1-thickness,0,h};


Circle(1) = { 1,2,3};
Circle(2) = { 3,2,1};
Line(3) = { 4,5};
Line(4) = { 5,6};
Line(5) = { 6,7};
Line(6) = { 7,4};
Line Loop(1) = {1,2};
Line Loop(2) = {3,4,5,6};
Plane Surface(3) = {2,1};

Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 6};
Line(10) = {5, 11};
Line(11) = {11, 10};
Line(12) = {10, 4};
Line Loop(13) = {8, 9, 5, 7};
Plane Surface(14) = {13};
Line Loop(15) = {3, 10, 11, 12};
Plane Surface(16) = {15};

Physical Line("rightTE") = {4};
Physical Line("interface") = {3,5};
Physical Line("leftTE") = {6};
Physical Line("cylinder") = {1,2};

Physical Line("leftHT") = {7,12};
Physical Line("rightHT") = {9,10};
Physical Line("outerHT") = {8,11};

Physical Surface("OmegaTE") = {3};
Physical Surface("OmegaHT") = {14,16};
