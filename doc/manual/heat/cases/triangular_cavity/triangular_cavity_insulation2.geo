h=0.025;

T = 0.00125;

//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {4, 0, 0, h};
//+
Point(3) = {4, 1, 0, h};
//+
Point(4) = {0, 1, 0, h};
//+
Point(5) = {5, 1, 0, h};
//+
Point(6) = {5, 4, 0, h};
//+
Point(7) = {4, 4, 0, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {3, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 3};
//+
Line(9) = {4, 7};
//+
Point(9) = {-0.6, 1.8, 0, h};
//+
Point(10) = {3.4, 4.8, 0, h};
//+
Line(10) = {9, 4};
//+
Line(11) = {9, 10};
//+
Line(12) = {10, 7};
//+
Point(11) = {T, 1, 0, h};
//+
Point(12) = {0, 1-T, 0, h};
//+
Point(13) = {T, 1-T, 0, h};
//+
Point(14) = {4-T, 1, 0, h};
//+
Point(15) = {4-T, 1-T, 0, h};
//+
Point(16) = {4, 1-T, 0, h};
//+
Point(17) = {4+T, 1, 0, h};
//+
Point(18) = {4+T, 1+T, 0, h};
//+
Point(19) = {4, 1+T, 0, h};
//+
Point(20) = {4, 4-T, 0, h};
//+
Point(21) = {4+T, 4-T, 0, h};
//+
Point(22) = {4+T, 4, 0, h};
//+
Point(24) = {4-0.6*T, 4+0.8*T, 0,h};
//+
Point(25) = {-0.6*T, 1+0.8*T, 0, h};
//+
Point(26) = {4*T, 1 +3*T, 0, h};
//+
Point(27) = {4-4*T, 4-3*T, 0, h};
//+
Point(28) = {3.4*T, 1+3.8*T, 0, h};
//+
Point(29) = {4-4.6*T, 4-2.2*T, 0, h};
//+
Split Curve {9} Point {26};
//+
Split Curve {10} Point {25};
//+
Split Curve {3} Point {11, 14};
//+
Split Curve {4} Point {12};
//+
Split Curve {2} Point {16};
//+
Split Curve {5} Point {17};
//+
Split Curve {8} Point {19, 20};
//+
Split Curve {7} Point {22};
//+
Split Curve {12} Point {24};
//+
Split Curve {13} Point {26, 27};
//+
Line(23) = {25, 28};
//+
Line(24) = {28, 26};
//+
Line(25) = {11, 13};
//+
Line(26) = {13, 12};
//+
Line(27) = {14, 15};
//+
Line(28) = {15, 16};
//+
Line(29) = {17, 18};
//+
Line(30) = {18, 19};
//+
Line(31) = {27, 29};
//+
Line(32) = {29, 24};
//+
Line(33) = {20, 21};
//+
Line(34) = {21, 22};
//+
Recursive Delete {
  Curve{22}; 
}
//+
Line(35) = {4, 26};
//+
Line(36) = {26, 27};
//+
Line(37) = {27, 7};
//+
Recursive Delete {
  Curve{14}; 
}
//+
Recursive Delete {
  Curve{21}; 
}
//+
Line(38) = {9, 25};
//+
Line(39) = {25, 4};
//+
Line(40) = {10, 24};
//+
Line(41) = {24, 7};
//+
Curve Loop(1) = {11, 40, -32, -31, -36, -24, -23, -38};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {39, 35, -24, -23};
//+
Curve Loop(3) = {31, 32, 41, -37};
//+
Plane Surface(2) = {2};
//+
Plane Surface(21) = {3};

//+
Recursive Delete {
  Curve{15}; 
}
//+
Recursive Delete {
  Curve{19}; 
}
//+
Recursive Delete {
  Curve{20}; 
}
//+
Recursive Delete {
  Curve{18}; 
}
//+
Recursive Delete {
  Curve{17}; 
}
//+
Recursive Delete {
  Curve{16}; 
}
//+
Line(42) = {4, 12};
//+
Line(43) = {12, 1};
//+
Line(44) = {4, 11};
//+
Line(45) = {11, 14};
//+
Line(46) = {2, 16};
//+
Line(47) = {17, 5};
//+
Line(48) = {7, 22};
//+
Line(49) = {22, 6};
//+
Line(50) = {7, 20};
//+
Line(51) = {20, 19};
//+
Point(30) = {4, 1, 0, h};
//+
Line(52) = {14, 30};
//+
Line(53) = {30, 17};
//+
Line(54) = {19, 30};
//+
Line(55) = {30, 16};
//+
Curve Loop(4) = {45, 27, 28, -46, -1, -43, -26, -25};
//+
Plane Surface(3) = {4};
//+
Curve Loop(5) = {52, 55, -28, -27};
//+
Curve Loop(6) = {44, 25, 26, -42};
//+
Plane Surface(4) = {5};
//+
Plane Surface(41) = {6};
//+
Curve Loop(7) = {29, 30, 54, 53};
//+
Curve Loop(8) = {34, -48, 50, 33};
//+
Plane Surface(5) = {7};
//+
Plane Surface(51) = {8};
//+
Curve Loop(9) = {49, -6, -47, 29, 30, -51, 33, 34};
//+
Plane Surface(6) = {9};
//+
Physical Curve("Fixed_T", 56) = {1};
//+
Physical Curve("FixedQ2", 57) = {6};
//+
Physical Curve("FixedQ3", 58) = {11};
//+
Physical Curve("RadiativeSurface1", 59) = {45};//, 44, 52};
//+
Physical Curve("RadiativeSurface2", 60) = {51}; //, 54, 50};
//+
Physical Curve("RadiativeSurface3", 61) = {36};//, 35, 37};
//+
Physical Surface("Rectangle1", 62) = {3};
//+
Physical Surface("Rectangle2", 63) = {6};
//+
Physical Surface("Rectangle3", 64) = {1};
//+
Physical Surface("Insulation1", 65) = {4};
//
Physical Surface("Insulation11", 651) = {41};
//+
Physical Surface("Insulation2", 66) = {5};
//+
Physical Surface("Insulation21", 661) = {51};
//+
Physical Surface("Insulation3", 67) = {2};

Physical Surface("Insulation31", 671) = {21};
