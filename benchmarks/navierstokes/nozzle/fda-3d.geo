hmax = 0.0008;
hmin = 0.00017;
lo = 0.18;  //0.18;0.144 

Mesh.CharacteristicLengthMin=hmin;
Mesh.CharacteristicLengthMax=hmax;
Mesh.CharacteristicLengthExtendFromBoundary=1;
Mesh.CharacteristicLengthFromPoints=1;

Point(1) = {-0.006, 0, -0.182685, hmax};
Point(2) = {-0.006, 0, -0.088, hmax};
Point(3) = {-0.006, 0, -0.064, hmax};
Point(4) = {-0.006, 0, -0.062685, hmax};
Point(5) = {-0.003410623759, 0, -0.048, hmin};
Point(6) = {-0.002, 0, -0.04, hmin};
Point(7) = {-0.002, 0, -0.02, hmin};
Point(8) = {-0.002, 0, -0.008, hmin};
Point(9) = {-0.002, 0, 0, hmin};
Point(10) = {-0.006, 0, 0, hmax};
Point(11) = {-0.006, 0, 0.008, hmax};
Point(12) = {-0.006, 0, 0.016, hmax};
Point(13) = {-0.006, 0, 0.024, hmax};
Point(14) = {-0.006, 0, 0.032, hmax};
Point(15) = {-0.006, 0, 0.06, hmax};
Point(16) = {-0.006, 0, 0.08, hmax};
Point(17) = {-0.006, 0, lo, hmax};


Point(18) = {0, 0, -0.182685, hmax};
Point(19) = {0, 0, -0.088, hmax};
Point(20) = {0, 0, -0.064, hmax};
Point(21) = {0, 0, -0.062685, hmax};
Point(22) = {0, 0, -0.048, hmax};
Point(23) = {0, 0, -0.04, hmin};
Point(24) = {0, 0, -0.02, hmin};
Point(25) = {0, 0, -0.008, hmin};
Point(26) = {0, 0, 0, hmin};
//Point(27) = {0, 0, 0, hmin};
Point(28) = {0, 0, 0.008, hmax};
Point(29) = {0, 0, 0.016, hmax};
Point(30) = {0, 0, 0.024, hmax};
Point(31) = {0, 0, 0.032, hmax};
Point(32) = {0, 0, 0.06, hmax};
Point(33) = {0, 0, 0.08, hmax};
Point(34) = {0, 0, lo, hmax};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 34};
Line(18) = {34, 33};
Line(19) = {33, 32};
Line(20) = {32, 31};
Line(21) = {31, 30};
Line(22) = {30, 29};
Line(23) = {29, 28};
Line(24) = {28, 26};
Line(26) = {26, 25};
Line(27) = {25, 24};
Line(28) = {24, 23};
Line(29) = {23, 22};
Line(30) = {22, 21};
Line(31) = {21, 20};
Line(32) = {20, 19};
Line(33) = {19, 18};
Line(34) = {18, 1};
Line(35) = {2, 19};
Line(36) = {3, 20};
Line(37) = {4, 21};
Line(38) = {5, 22};
Line(39) = {6, 23};
Line(40) = {7, 24};
Line(41) = {8, 25};
Line(42) = {9, 26};
Line(43) = {11, 28};
Line(44) = {12, 29};
Line(45) = {13, 30};
Line(46) = {14, 31};
Line(47) = {15, 32};
Line(48) = {16, 33};

Line Loop(49) = {34, 1, 35, 33};
Plane Surface(50) = {49};
Line Loop(51) = {35, -32, -36, -2};
Plane Surface(52) = {51};
Line Loop(53) = {36, -31, -37, -3};
Plane Surface(54) = {53};
Line Loop(55) = {37, -30, -38, -4};
Plane Surface(56) = {55};
Line Loop(57) = {38, -29, -39, -5};
Plane Surface(58) = {57};
Line Loop(59) = {39, -28, -40, -6};
Plane Surface(60) = {59};
Line Loop(61) = {40, -27, -41, -7};
Plane Surface(62) = {61};
Line Loop(63) = {41, -26, -42, -8};
Plane Surface(64) = {63};
Line Loop(65) = {42, -24, -43, -10, -9};
Plane Surface(66) = {65};
Line Loop(67) = {43, -23, -44, -11};
Plane Surface(68) = {67};
Line Loop(69) = {44, -22, -45, -12};
Plane Surface(70) = {69};
Line Loop(71) = {45, -21, -46, -13};
Plane Surface(72) = {71};
Line Loop(73) = {46, -20, -47, -14};
Plane Surface(74) = {73};
Line Loop(75) = {47, -19, -48, -15};
Plane Surface(76) = {75};
Line Loop(77) = {48, -18, -17, -16};
Plane Surface(78) = {77};

Physical Surface(79) = {50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{96, 113, 130, 147, 164, 181, 198, 215, 237, 254, 271, 288, 305, 322, 339};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{356, 373, 390, 407, 424, 441, 458, 475, 497, 514, 531, 548, 565, 582, 599};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{616, 633, 650, 667, 684, 701, 718, 735, 757, 774, 791, 808, 825, 842, 859};
}
Physical Surface("inlet") = {607, 347, 87, 867};
Physical Surface("face1") = {614, 354, 94, 874};
Physical Surface("face2") = {628, 368, 108, 886};
Physical Surface("face3") = {918, 662, 402, 142};
Physical Surface("face4") = {696, 436, 176, 950};
Physical Surface("face5") = {713, 453, 193, 966};
Physical Surface("face6") = {982, 730, 470, 210};
Physical Surface("face7") = {748, 999, 228, 488};
Physical Surface("face8") = {1019, 249, 509, 769};
Physical Surface("face9") = {786, 526, 266, 1035};
Physical Surface("face10") = {803, 543, 283, 1051};
Physical Surface("face11") = {820, 560, 300, 1067};
Physical Surface("face12") = {837, 577, 317, 1083};
Physical Surface("outlet") = {854, 594, 334, 1099};
Physical Surface("wall") = {871, 611, 351, 91, 632, 372, 112, 890, 649, 389, 129, 906, 922, 666, 406, 146, 938, 683, 423, 163, 954, 700, 440, 180, 970, 717, 457, 197, 986, 734, 474, 214, 756, 1007, 236, 496, 1003, 752, 492, 232, 1023, 773, 513, 253, 1039, 790, 530, 270, 1055, 807, 547, 287, 824, 564, 304, 1071, 1087, 581, 841, 321, 1103, 858, 598, 338};
Physical Volume("fluid") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60};
