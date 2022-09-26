h=0.3;

// Face arriere

Point(1) = {-2,1,0,h};
Point(2) = {-2,-1,0,h};
Point(3) = {0,-1,0,h};
Point(4) = {0,1,0,h};
Point(5) = {2,1,0,h};
Point(6) = {2,-1,0,h};

Line(1) = {4,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,3};

Line Loop(21) = {1,2,3,4};
Line Loop(22) = {-5,-4,-7,-6};

Plane Surface(81) = {21};
Plane Surface(82) = {22};

// Face avant

Point(7) = {-2,1,2,h};
Point(8) = {-2,-1,2,h};
Point(9) = {0,-1,2,h};
Point(10) = {0,1,2,h};
Point(11) = {2,1,2,h};
Point(12) = {2,-1,2,h};

Line(8) = {10,7};
Line(9) = {7,8};
Line(10) = {8,9};
Line(11) = {9,10};
Line(12) = {10,11};
Line(13) = {11,12};
Line(14) = {12,9};

Line Loop(23) = {8,9,10,11};
Line Loop(24) = {-12,-11,-14,-13};

Plane Surface(83) = {23};
Plane Surface(84) = {24};

// Profondeur

Line(15) = {1,7};
Line(16) = {2,8};
Line(17) = {4,10};
Line(18) = {3,9};
Line(19) = {5,11};
Line(20) = {6,12};

Line Loop(25) = {15,-8,-17,1};
Line Loop(26) = {15,9,-16,-2};
Line Loop(27) = {16,10,-18,-3};

Plane Surface(85) = {25};
Plane Surface(86) = {26};
Plane Surface(87) = {27};

Line Loop(28) = {17,12,-19,-5};
Line Loop(29) = {19,13,-20,-6};
Line Loop(30) = {18,-14,-20,7};

Plane Surface(88) = {28};
Plane Surface(89) = {29};
Plane Surface(90) = {30};

Line Loop(31) = {17,-11,-18,4};
Plane Surface(91) = {31};

// Cube 1

Surface Loop(1) = {85,86,87,91,81,83};
Volume(1) = {1};

// Cube 2

Surface Loop(2) = {88,91,90,89,82,84};
Volume(2) = {2};

Physical Surface("Gamma_1",101) = {85,81,83,87};
Physical Surface("Gamma_2",102) = {88,90,82,84};
Physical Surface("Gamma_dir",103) = {86};
Physical Surface("Gamma_neu",104) = {89};
Physical Surface("Gamma_int",105) = {91};

Physical Volume("Omega_1",201) = {1};
Physical Volume("Omega_2",202) = {2};
