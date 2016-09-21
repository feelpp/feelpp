cl__1 = 1e+22;

h = 0.5;

Mesh.CharacteristicLengthFactor=h;
Mesh.CharacteristicLengthMax=0.8;


Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(5) = {-1, 0, 0, 1.0};
Point(6) = {0, -1, 0, 1.0};
Circle(1) = {6, 1, 2};
Circle(2) = {2, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 6};
Point(7) = {0.4, 0, 0, 1.0};
Point(8) = {-0.4, 0, 0, 1.0};
Point(9) = {0, 0.4, 0, 1.0};
Point(10) = {0, -0.4, 0, 1.0};
Circle(5) = {8, 1, 9};
Circle(6) = {9, 1, 7};
Circle(7) = {7, 1, 10};
Circle(8) = {10, 1, 8};
Line Loop(9) = {3, 4, 1, 2};
Line Loop(10) = {5, 6, 7, 8};
Plane Surface(11) = {9, 10};
Plane Surface(12) = {10};


Plane Surface(13) = {10};
Extrude {0, 0, 3} {
  Surface{12};
}

Line Loop(36) = {20, -18, -29, 8};
Ruled Surface(37) = {36};
Line Loop(38) = {7, 29, -17, -25};
Ruled Surface(39) = {38};
Line Loop(40) = {16, -25, -6, 21};
Ruled Surface(41) = {40};
Line Loop(42) = {15, -21, -5, 20};
Ruled Surface(43) = {42};
Line Loop(44) = {15, 16, 17, 18};
Plane Surface(45) = {44};
Physical Surface("lateral_internal") = {34, 30, 26, 22};
Physical Surface("top_internal") = {35};
Physical Surface("bottom_internal") = {12};
Extrude {0, 0, 3} {
  Surface{11};
}
Line Loop(88) = {57, -47, -56, 3};
Ruled Surface(89) = {88};
Line Loop(90) = {50, -56, -2, 65};
Ruled Surface(91) = {90};
Line Loop(92) = {49, -65, -1, 61};
Ruled Surface(93) = {92};
Line Loop(94) = {48, -61, -4, 57};
Ruled Surface(95) = {94};
Line Loop(96) = {47, 48, 49, 50};
Plane Surface(97) = {44, 96};
Physical Surface("lateral_external") = {62, 58, 70, 66};
Physical Surface("top_external") = {87};
Physical Surface("bottom_external") = {11};
// Surface Loop(98) = {11, 87, 70, 62, 58, 66, 22, 26, 30, 34};
// Volume(100) = {98};
Physical Volume("external_cylinder") = {2};
// Surface Loop(101) = {35, 22, 26, 30, 34, 12};
// Volume(102) = {101};
Physical Volume("internal_cylinder") = {1};

