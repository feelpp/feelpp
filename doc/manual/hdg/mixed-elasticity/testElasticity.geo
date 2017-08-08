cl__1 = 1e+22;

h=0.1;

Mesh.CharacteristicLengthFactor=h;
Mesh.CharacteristicLengthMax=0.8;

Point(1) = {0.095, 0, -0.015, 1e+22};
Point(2) = {0, 0, -0.015, 1e+22};
Point(3) = {0.04750000000000001, 0.08227241335952167, -0.015, 1e+22};
Point(4) = {-0.08227241335952169, -0.04749999999999997, -0.015, 1e+22};
Point(5) = {0.008, 0, -0.015, 1e+22};
Point(6) = {0, 0.008, -0.015, 1e+22};
Point(7) = {-0.008, 0, -0.015, 1e+22};
Point(8) = {0, -0.008, -0.015, 1e+22};
Point(9) = {0.095, 0, 0.015, 1e+22};
Point(10) = {0, 0, 0.015, 1e+22};
Point(11) = {0.04750000000000001, 0.08227241335952167, 0.015, 1e+22};
Point(14) = {-0.08227241335952169, -0.04749999999999997, 0.015, 1e+22};
Point(15) = {0, 0.008, 0.015, 1e+22};
Point(17) = {0.008, 0, 0.015, 1e+22};
Point(20) = {0, -0.008, 0.015, 1e+22};
Point(23) = {-0.008, 0, 0.015, 1e+22};
Circle(1) = {1, 2, 3};
Circle(2) = {3, 2, 4};
Circle(3) = {4, 2, 1};
Circle(5) = {5, 2, 6};
Circle(6) = {6, 2, 7};
Circle(7) = {7, 2, 8};
Circle(8) = {8, 2, 5};
Circle(9) = {9, 10, 11};
Line(10) = {1, 9};
Line(11) = {3, 11};
Circle(13) = {11, 10, 14};
Line(15) = {4, 14};
Circle(17) = {14, 10, 9};
Circle(21) = {15, 10, 17};
Line(22) = {6, 15};
Line(23) = {5, 17};
Circle(25) = {17, 10, 20};
Line(27) = {8, 20};
Circle(29) = {20, 10, 23};
Line(31) = {7, 23};
Circle(33) = {23, 10, 15};
Line Loop(12) = {1, 11, -9, -10};
Ruled Surface(12) = {12};
Line Loop(16) = {2, 15, -13, -11};
Ruled Surface(16) = {16};
Line Loop(20) = {3, 10, -17, -15};
Ruled Surface(20) = {20};
Line Loop(24) = {-5, 23, -21, -22};
Ruled Surface(24) = {24};
Line Loop(28) = {-8, 27, -25, -23};
Ruled Surface(28) = {28};
Line Loop(32) = {-7, 31, -29, -27};
Ruled Surface(32) = {32};
Line Loop(36) = {-6, 22, -33, -31};
Ruled Surface(36) = {36};
Line Loop(1000) = {1, 2, 3, 6, 7, 8, 5};
Plane Surface(1000) = {1000};
Line Loop(1001) = {9, 13, 17, -33, -29, -25, -21};
Plane Surface(1001) = {1001};
Surface Loop(3000) = {12, 16, 20, 24, 28, 32, 36, 1000, 1001};
Volume(3000) = {3000};

// Physical labels
Physical Surface("out") = {12, 16, 20}; // Integral
// Physical Surface("top") = {1000} ; // Neumann
// Physical Surface("bottom") = {1001} ; // Neumann
// Physical Surface("in") = {24, 28, 32, 36}; // Dirichlet
// Physical Volume("omega") = {3000};

Delete {
  Surface{36};
}
Delete {
  Surface{32};
}
Delete {
  Point{10};
}
Delete {
  Line{21};
}
Delete {
  Line{25, 21, 29, 33, 31, 6, 5, 8, 23, 7};
}
Delete {
  Volume{3000};
}
Delete {
  Surface{36, 24, 28, 32};
}
Delete {
  Line{29, 21, 22, 33, 31, 7, 6, 5, 27, 23, 8, 25};
}
Delete {
  Line{29};
}
Delete {
  Line{7, 8, 5, 6};
}
Delete {
  Surface{1001};
}
Delete {
  Surface{1000};
}
Delete {
  Line{33, 29};
}
Delete {
  Line{25};
}
Delete {
  Line{21};
}
Delete {
  Line{5};
}
Delete {
  Line{6};
}
Delete {
  Line{7};
}
Delete {
  Line{8};
}
Delete {
  Point{10, 20, 17, 15};
}
Delete {
  Point{6};
}
Delete {
  Point{10};
}
Delete {
  Point{2};
}
Delete {
  Point{8, 7, 2, 5};
}
Delete {
  Point{23};
}
Delete {
  Point{10};
}
Point(23) = {0.1, 0.1, 0, 1.0};
Delete {
  Point{23};
}
Point(23) = {0.001, 0.001, 0, 1.0};
Delete {
  Point{23};
}
Point(23) = {0, 0, -0, 1.0};
Point(24) = {0.02, 0.02, 0, 1.0};
Point(25) = {0.02, 0.02, -0.015, 1.0};
Point(26) = {0.02, 0.02, 0.015, 1.0};

Point(27) = {0.04, 0.04, -0.015, 1.0};
Point(28) = {0.04, 0.04, 0.015, 1.0};
Line(3001) = {10, 2};
Line(3002) = {28, 27};
Circle(3003) = {10, 26, 28};
Circle(3004) = {28, 26, 10};
Line Loop(3005) = {3003, 3004};

Circle(3006) = {2, 25, 27};
Circle(3007) = {27, 25, 2};
Line Loop(3008) = {17, 9, 13, 17, 9};

Line Loop(3009) = {3003, 3002, -3006, -3001};
Ruled Surface(3010) = {3009};
Line Loop(3011) = {3002, 3007, -3001, -3004};
Ruled Surface(3012) = {3011};
Line Loop(3013) = {17, 9, 13};
Plane Surface(3014) = {3005, 3013};
Line Loop(3015) = {3, 1, 2};
Line Loop(3016) = {3006, 3007};
Plane Surface(3017) = {3015, 3016};

// Physical Surface("top") = {1000} ; // Neumann
// Physical Surface("bottom") = {1001} ; // Neumann
// Physical Surface("in") = {24, 28, 32, 36}; // Dirichlet

Physical Surface("top") = {3017}; 	// Neumann
Physical Surface("bottom") = {3014};	// Neumann
Physical Surface("in") = {3010, 3012};	// Dirichlet


Surface Loop(3018) = {3014, 3010, 3012, 3017, 20, 12, 16};
Volume(3019) = {3018};
Physical Volume("omega") = {3019};
