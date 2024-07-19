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


Point(12) = {1+0.5,1+thickness,0,h};
Point(13) = {1+0.5,-1-thickness,0,h};
Point(14) = {1+1,1+thickness,0,h};
Point(15) = {1+1,-1-thickness,0,h};
Point(16) = {1+1+1,1+thickness,0,h};
Point(17) = {1+1+1,-1-thickness,0,h};



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

Line(13) = {9, 12};
Line(14) = {12, 14};
Line(15) = {14, 16};
Line(16) = {16, 17};
Line(17) = {17, 15};
Line(18) = {15, 13};
Line(19) = {13, 11};
Line(20) = {12, 13};
Line(21) = {15, 14};
Curve Loop(16) = {9, -4, 10, -19, -20, -13};
Plane Surface(17) = {16};
Curve Loop(17) = {14, -21, 18, -20};
Plane Surface(18) = {17};
Curve Loop(18) = {21, 15, 16, 17};
Plane Surface(19) = {18};



Physical Line("rightTE") = {4};
Physical Line("interface") = {3,5};
Physical Line("leftTE") = {6};
Physical Line("cylinder") = {1,2};

Physical Line("leftHT") = {7,12};
Physical Line("rightHT") = {9,10};
Physical Line("outerHT") = {8,11};

Physical Line("topTE2") = {15};
Physical Line("bottomTE2") = {17};
Physical Line("rightTE2") = {16};


Physical Surface("OmegaTE") = {3};
Physical Surface("OmegaHT") = {14,16};

Physical Surface("OmegaHT2") = {17};
Physical Surface("OmegaHT3") = {18};
Physical Surface("OmegaTE2") = {19};


