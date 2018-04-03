h=0.001;

W=0.5;
w=0.015;
h1=0.006;
h2=0.005;
h3=0.0415;
h4=0.0365;
t4=0.0015;

Point(1) = {0,0,0.0,h};
Point(2) = {W,0,0.0,h};
Point(3) = {W,h3+h1,0.0,h};
Point(4) = {0,h3+h1,0.0,h};

Point(5) = {W,t4,0.0,h};
Point(6) = {t4,t4,0.0,h};
Point(7) = {t4,h3-(h2+t4),0.0,h};
Point(8) = {w,h3-(h2+t4),0.0,h};
Point(9) = {w,h3-(h2),0.0,h};
Point(10) = {0,h3-(h2),0.0,h};

Point(11) = {0,h3,0.0,h};
Point(12) = {w,h3,0.0,h};
Point(13) = {W,h3,0.0,h};
Line(1) = {2, 1};
Line(2) = {1, 10};
Line(3) = {10, 11};
Line(4) = {11, 4};
Line(5) = {4, 3};
Line(6) = {3, 13};
Line(7) = {13, 5};
Line(8) = {5, 2};
Line(9) = {5, 6};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 9};
Line(13) = {9, 10};
Line(14) = {9, 12};
Line(15) = {12, 11};
Line(16) = {12, 13};
Line Loop(17) = {5, 6, -16, 15, 4};
Plane Surface(18) = {17};
Line Loop(19) = {15, -3, -13, 14};
Plane Surface(20) = {19};
Line Loop(21) = {16, 7, 9, 10, 11, 12, 14};
Plane Surface(22) = {21};
Line Loop(23) = {9, 10, 11, 12, 13, -2, -1, -8};
Plane Surface(24) = {23};

Physical Line("top") = {5};
Physical Line("bottom") = {1};
Physical Surface("mat1") = {18};
Physical Surface("mat2") = {20};
Physical Surface("mat3") = {22};
Physical Surface("mat4") = {24};
