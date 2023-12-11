h=0.0125;

T = 0.0025;

L = 0.025;

//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 3, L, 0};
//+
Rectangle(2) = {3, L, -0, L, 3, 0};
//+
Rectangle(3) = {-L, L, -0, L, 3, 0};
//+
Rectangle(4) = {0, 3+L, 0, 3, L, 0};
//+
Point(17) = {0+T, L, 0, h};
//+
Point(18) = {0, L-T, 0, h};
//+
Point(19) = {0, L+T, 0, h};
//+
Point(20) = {-T, L, 0, h};
//+
Point(21) = {T, L-T, 0, h};
//+
Point(22) = {-T, L+T, 0, h};
//+
Point(23) = {-T, 3+L, 0, h};
//+
Point(24) = {-0, 3+L-T, 0, h};
//+
Point(25) = {-T, 3+L-T, 0, h};
//+
Point(26) = {0, 3+L+T, 0, h};
//+
Point(27) = {T, 3+L+T, 0, h};
//+
Point(28) = {T, 3+L, 0, h};
//+
Point(29) = {3-T, 3+L, 0, h};
//+
Point(30) = {3-T, 3+L+T, 0, h};
//+
Point(31) = {3, 3+L+T, 0, h};
//+
Point(32) = {3+T, 3+L, 0, h};
//+
Point(33) = {3, 3+L-T, 0, h};
//+
Point(34) = {3+T, 3+L-T, 0, h};
//+
Point(35) = {3, L+T, 0, h};
//+
Point(36) = {3+T, L+T, 0, h};
//+
Point(37) = {3+T, L, 0, h};
//+
Recursive Delete {
  Curve{3}; 
}
//+
Recursive Delete {
  Curve{3}; 
}
//+
Delete {
  Surface{4}; Surface{3}; Surface{2}; Surface{1}; 
}
//+
Point(38) = {3-T, L, 0, h};
//+
Point(39) = {3-T, L-T, 0, h};
//+
Point(40) = {3, L-T, 0, h};
//+
Delete {
  Curve{10}; Curve{11}; Curve{9}; Curve{4}; Curve{3}; Curve{2}; Curve{5}; Curve{8}; Curve{7}; Curve{14}; Curve{13}; Curve{16}; 
}
//+
Line(16) = {12, 23};
//+
Line(17) = {23, 11};
//+
Line(18) = {11, 26};
//+
Line(19) = {26, 16};
//+
Line(20) = {23, 25};
//+
Line(21) = {25, 24};
//+
Line(22) = {24, 11};
//+
Line(23) = {11, 28};
//+
Line(24) = {27, 28};
//+
Line(25) = {26, 27};
//+
Line(26) = {28, 29};
//+
Line(27) = {29, 30};
//+
Line(28) = {30, 31};
//+
Line(29) = {31, 15};
//+
Line(30) = {31, 8};
//+
Line(31) = {29, 8};
//+
Line(32) = {8, 32};
//+
Line(33) = {32, 7};
//+
Line(34) = {8, 33};
//+
Line(35) = {33, 34};
//+
Line(36) = {34, 32};
//+
Line(37) = {33, 35};
//+
Line(38) = {35, 36};
//+
Line(39) = {36, 37};
//+
Line(40) = {37, 3};
//+
Line(41) = {3, 35};
//+
Line(42) = {37, 6};
//+
Line(43) = {38, 3};
//+
Line(44) = {3, 40};
//+
Line(45) = {40, 2};
//+
Line(46) = {38, 39};
//+
Line(47) = {39, 40};
//+
Line(48) = {4, 17};
//+
Line(49) = {17, 38};
//+
Line(50) = {4, 18};
//+
Line(51) = {18, 21};
//+
Line(52) = {21, 17};
//+
Line(53) = {18, 1};
//+
Line(54) = {9, 20};
//+
Line(55) = {20, 4};
//+
Line(56) = {4, 19};
//+
Line(57) = {19, 22};
//+
Line(58) = {22, 20};
//+
Line(59) = {19, 24};
//+
Curve Loop(5) = {15, -19, 25, 24, 26, 27, 28, 29};
//+
Plane Surface(1) = {5};
//+
Curve Loop(6) = {27, 28, 30, -31};
//+
Plane Surface(2) = {6};
//+
Curve Loop(7) = {32, -36, -35, -34};
//+
Plane Surface(3) = {7};
//+
Curve Loop(8) = {33, -6, -42, -39, -38, -37, 35, 36};
//+
Plane Surface(4) = {8};
//+
Curve Loop(9) = {38, 39, 40, 41};
//+
Plane Surface(5) = {9};
//+
Curve Loop(10) = {43, 44, -47, -46};
//+
Plane Surface(6) = {10};
//+
Curve Loop(11) = {45, -1, -53, 51, 52, 49, 46, 47};
//+
Plane Surface(7) = {11};
//+
Curve Loop(12) = {51, 52, -48, 50};
//+
Plane Surface(8) = {12};
//+
Curve Loop(13) = {55, 56, 57, 58};
//+
Plane Surface(9) = {13};
//+
Curve Loop(14) = {59, -21, -20, -16, 12, 54, -58, -57};
//+
Plane Surface(10) = {14};
//+
Curve Loop(15) = {21, 22, -17, 20};
//+
Plane Surface(11) = {15};
//+
Curve Loop(16) = {23, -24, -25, -18};
//+
Plane Surface(12) = {16};
//+
Physical Surface("Rectangle1", 60) = {7};
//+
Physical Surface("Rectangle2", 61) = {4};
//+
Physical Surface("Rectangle3", 62) = {1};
//+
Physical Surface("Rectangle4", 63) = {10};

Characteristic Length{ PointsOf{ Surface{:}; } } = h;
//+
Physical Curve("ExternalSurface4", 64) = {12};
//+
Physical Curve("ExternalSurface3", 65) = {15};
//+
Physical Curve("ExternalSurface2", 66) = {6};
//+
Physical Curve("ExternalSurface1", 67) = {1};
//+
Physical Curve("RadiativeSurface1", 68) = {49};
//+
Physical Curve("RadiativeSurface2", 69) = {37};
//+
Physical Curve("RadiativeSurface3", 70) = {26};
//+
Physical Curve("RadiativeSurface4", 71) = {59};
//+
Physical Surface("Insulation1", 72) = {8};
//+
Physical Surface("Insulation2", 73) = {6};
//+
Physical Surface("Insulation3", 74) = {5};
//+
Physical Surface("Insulation4", 75) = {3};
//+
Physical Surface("Insulation5", 76) = {2};
//+
Physical Surface("Insulation6", 77) = {12};
//+
Physical Surface("Insulation7", 78) = {11};
//+
Physical Surface("Insulation8", 79) = {9};
