h=0.1;
W=1;
w=0.4;
Point(1) = {0,0,0.0,h};
Point(2) = {W,0,0.0,h};
Point(3) = {0,w,0.0,h};
Point(4) = {W,w,0.0,h};
Point(5) = {0,2*w,0.0,h};
Point(6) = {W,2*w,0.0,h};
Point(7) = {0,3*w,0.0,h};
Point(8) = {W,3*w,0.0,h};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line(5) = {4, 6};
Line(6) = {6, 5};
Line(7) = {5, 3};
Line(8) = {6, 8};
Line(9) = {8, 7};
Line(10) = {7, 5};
Line Loop(11) = {1, 2, 3, 4};
Plane Surface(12) = {11};
Line Loop(13) = {3, -7, -6, -5};
Plane Surface(14) = {13};
Line Loop(15) = {6, -10, -9, -8};
Plane Surface(16) = {15};

Physical Surface("mat1") = {12};
Physical Surface("mat2") = {14};
Physical Surface("mat3") = {16};

Physical Line("line-h1") = {1};
Physical Line("line-h2") = {3};
Physical Line("line-h3") = {6};
Physical Line("line-h4") = {9};

Physical Line("line-vl1") = {4};
Physical Line("line-vl2") = {7};
Physical Line("line-vl3") = {10};

Physical Line("line-vr1") = {2};
Physical Line("line-vr2") = {5};
Physical Line("line-vr3") = {8};
