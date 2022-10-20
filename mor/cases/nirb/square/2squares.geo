h=0.1;

Point(1) = {-2,1,0,h};
Point(2) = {-2,-1,0,h};
Point(3) = {0,-1,0,h};
Point(4) = {0,1,0,h};
Point(5) = {2,1,0,h};
Point(6) = {2,-1,0,h};

Line(11) = {4,1};
Line(12) = {1,2};
Line(13) = {2,3};
Line(14) = {3,4};
Line(15) = {4,5};
Line(16) = {5,6};
Line(17) = {6,3};

Line Loop(21) = {11,12,13,14};
Line Loop(22) = {-15,-14,-17,-16};

Plane Surface(81) = {21};
Plane Surface(82) = {22};

Physical Line("Gamma_1",91) = {11,13};
Physical Line("Gamma_2",92) = {15,17};
Physical Line("Gamma_dir",93) = {12};
Physical Line("Gamma_int",94) = {14};
Physical Line("Gamma_neu",95) = {16};

Physical Surface("Omega_1",101) = {81};
Physical Surface("Omega_2",102) = {82};