h=0.025;

//floor
Point(1) = {0, 0, 0, h};
Point(2) = {1.5, 0, 0, h};
Point(3) = {1.5, 0, 0.5, h};
Point(4) = {0, 0, 0.5, h};
//roof
Point(5) = {0, 1, 0, h};
Point(6) = {1.5, 1, 0, h};
Point(7) = {1.5, 1, 0.5, h};
Point(8) = {0, 1, 0.5, h};
//outlet
Point(9)  = {1.5, 0, 0.15, h};
Point(10) = {1.5, 0, 0.35, h};
Point(11) = {1.5, 0.10, 0.35, h};
Point(12) = {1.5, 0.10, 0.15, h};
//inlet
Point(13) = {0, 1, 0.15, h};
Point(14) = {0, 1, 0.35, h};
Point(15) = {0, 0.90, 0.35, h};
Point(16) = {0, 0.90, 0.15, h};
//cpu1
Point(17) = {0.40, 0, 0.15, h};
Point(18) = {0.50, 0, 0.15, h};
Point(19) = {0.50, 0, 0.35, h};
Point(20) = {0.40, 0, 0.35, h};
Point(21) = {0.40, 0.25, 0.15, h};
Point(22) = {0.50, 0.25, 0.15, h};
Point(23) = {0.50, 0.25, 0.35, h};
Point(24) = {0.40, 0.25, 0.35, h};

//cpu2
Point(25) = {0.70, 0, 0.15, h};
Point(26) = {0.90, 0, 0.15, h};
Point(27) = {0.90, 0, 0.35, h};
Point(28) = {0.70, 0, 0.35, h};
Point(29) = {0.70, 0.5, 0.15, h};
Point(30) = {0.90, 0.5, 0.15, h};
Point(31) = {0.90, 0.5, 0.35, h};
Point(32) = {0.70, 0.5, 0.35, h};

Line(1) = {4, 3};
Line(2) = {3, 10};
Line(3) = {10, 9};
Line(4) = {9, 2};
Line(5) = {2, 1};
Line(6) = {1, 4};
Line(7) = {8, 14};
Line(8) = {14, 13};
Line(9) = {13, 5};
Line(10) = {5, 6};
Line(11) = {6, 7};
Line(12) = {7, 8};
Line(13) = {8, 4};
Line(14) = {3, 7};
Line(15) = {6, 2};
Line(16) = {5, 1};
Line(17) = {10, 11};
Line(19) = {11, 12};
Line(20) = {12, 9};
Line(21) = {13, 16};
Line(22) = {16, 15};
Line(23) = {15, 14};
//cpu1
Line(24) = {20, 19};
Line(25) = {19, 18};
Line(26) = {18, 17};
Line(27) = {17, 20};
Line(28) = {24, 23};
Line(29) = {23, 22};
Line(30) = {22, 21};
Line(31) = {21, 24};
Line(32) = {20, 24};
Line(33) = {19, 23};
Line(34) = {18, 22};
Line(35) = {17, 21};
//cpu2
Line(36) = {28, 27};
Line(37) = {27, 26};
Line(38) = {26, 25};
Line(39) = {25, 28};
Line(40) = {32, 31};
Line(41) = {31, 30};
Line(42) = {30, 29};
Line(43) = {29, 32};
Line(44) = {28, 32};
Line(45) = {27, 31};
Line(46) = {26, 30};
Line(47) = {25, 29};

Line Loop(48) = {24, 33, -28, -32};
Plane Surface(49) = {48};
Line Loop(50) = {25, 34, -29, -33};
Plane Surface(51) = {50};
Line Loop(52) = {26, 35, -30, -34};
Plane Surface(53) = {52};
Line Loop(54) = {27, 32, -31, -35};
Plane Surface(55) = {54};
Line Loop(56) = {28, 29, 30, 31};
Plane Surface(57) = {56};
Line Loop(58) = {36, 45, -40, -44};
Plane Surface(59) = {58};
Line Loop(60) = {37, 46, -41, -45};
Plane Surface(61) = {60};
Line Loop(62) = {38, 47, -42, -46};
Plane Surface(63) = {62};
Line Loop(64) = {39, 44, -43, -47};
Plane Surface(65) = {64};
Line Loop(66) = {40, 41, 42, 43};
Plane Surface(67) = {66};
Line Loop(68) = {1, 2, 3, 4, 5, 6};
Line Loop(69) = {36, 37, 38, 39};
Line Loop(70) = {24, 25, 26, 27};
Plane Surface(71) = {68, 69, 70};
Line Loop(72) = {3, -20, -19, -17};
Plane Surface(73) = {72};
Line Loop(74) = {2, 17, 19, 20, 4, -15, 11, -14};
Plane Surface(75) = {74};
Line Loop(76) = {11, 12, 7, 8, 9, 10};
Plane Surface(77) = {76};
Line Loop(78) = {23, 8, 21, 22};
Plane Surface(79) = {78};
Line Loop(80) = {9, 16, 6, -13, 7, -23, -22, -21};
Plane Surface(81) = {80};
Line Loop(82) = {1, 14, 12, 13};
Plane Surface(83) = {82};
Line Loop(84) = {5, -16, 10, 15};
Plane Surface(85) = {84};

Surface Loop(86) = {81, 77, 75, 71, 83, 73, 85, 65, 59, 61, 63, 67, 55, 49, 51, 53, 57, 79};
Volume(87) = {86};

Physical Surface("inlet") = {79};
Physical Surface("outlet") = {73};
Physical Surface("cpu1") = {49, 51, 53, 55, 57};
Physical Surface("cpu2") = {59, 61, 63, 65, 67};
Physical Surface("wall") = {71, 75, 77, 81, 85, 83};
Physical Volume("omega") = {87};
