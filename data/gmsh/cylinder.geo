hmax = 0.0008;
hmin = 0.00017;
Point(5) = {0.006, 0, 0.18, hmax};
Point(6) = {0, 0, 0.18, hmax};
Point(7) = {0.006, 0, 0, hmax};
Point(8) = {0, 0, 0, hmin};
Point(10) = {0, 0, 0, hmin};
Point(27) = {0.006, 0, 0.008, hmax};
Point(28) = {0, 0, 0.008, hmax};
Point(29) = {0.006, 0, 0.016, hmax};
Point(30) = {0, 0, 0.016, hmax};
Point(31) = {0.006, 0, 0.024, hmax};
Point(32) = {0, 0, 0.024, hmax};
Point(33) = {0.006, 0, 0.032, hmax};
Point(34) = {0, 0, 0.032, hmax};
Point(35) = {0.006, 0, 0.06, hmax};
Point(36) = {0, 0, 0.06, hmax};
Point(37) = {0.006, 0, 0.08, hmax};
Point(38) = {0, 0, 0.08, hmax};
Line(1) = {8, 7};
Line(2) = {7, 27};
Line(3) = {27, 29};
Line(4) = {29, 31};
Line(5) = {31, 33};
Line(6) = {33, 35};
Line(7) = {35, 37};
Line(8) = {37, 5};
Line(9) = {5, 6};
Line(10) = {6, 38};
Line(11) = {38, 36};
Line(12) = {36, 34};
Line(13) = {34, 32};
Line(14) = {32, 30};
Line(15) = {30, 28};
Line(16) = {28, 8};
Line(17) = {28, 27};
Line(18) = {30, 29};
Line(19) = {32, 31};
Line(20) = {34, 33};
Line(21) = {36, 35};
Line(22) = {38, 37};
Line Loop(25) = {1, 2, -17, 16};
Plane Surface(25) = {25};
Line Loop(27) = {17, 3, -18, 15};
Plane Surface(27) = {27};
Line Loop(29) = {18, 4, -19, 14};
Plane Surface(29) = {29};
Line Loop(31) = {19, 5, -20, 13};
Plane Surface(31) = {31};
Line Loop(33) = {20, 6, -21, 12};
Plane Surface(33) = {33};
Line Loop(35) = {21, 7, -22, 11};
Plane Surface(35) = {35};
Line Loop(37) = {22, 8, 9, 10};
Plane Surface(37) = {37};
Physical Line(38) = {1};
Physical Line(39) = {2, 3, 4, 5, 6, 7, 8};
Physical Line(40) = {10, 11, 12, 13, 14, 15, 16};
Physical Line(41) = {9};
Physical Line(42) = {17};
Physical Line(43) = {18};
Physical Line(44) = {19};
Physical Line(45) = {20};
Physical Line(46) = {21};
Physical Line(47) = {22};
Physical Surface(48) = {25, 27, 29, 31, 33, 35, 37};


Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{25, 27, 29, 31, 33, 35, 37};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{65, 82, 99, 116, 133, 150, 167};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{184, 201, 218, 235, 252, 269, 286};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{303, 320, 337, 354, 371, 388, 405};
}
Physical Surface("inlet") = {56, 175, 294, 413};
Physical Surface("outlet") = {165, 284, 403, 510};
Physical Surface("face1") = {63, 182, 301, 420};
Physical Surface("face2") = {80, 199, 318, 435};
Physical Surface("face3") = {97, 216, 335, 450};
Physical Surface("face4") = {114, 233, 352, 465};
Physical Surface("face5") = {131, 250, 369, 480};
Physical Surface("face6") = {148, 267, 386, 495};
Physical Surface("wall") = {60, 77, 94, 111, 128, 162, 145, 179, 196, 213, 230, 247, 264, 281, 507, 492, 477, 462, 447, 432, 417, 400, 383, 366, 349, 332, 315, 298}; 
Physical Volume("fluid") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};

