h=0.10;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1.85, 0, 0, h};
Point(4) = {2.15, 0, 0, h};
Point(5) = {3, 0, 0, h};
Point(6) = {4, 0, 0, h};
Point(7) = {1.85, 0.3, 0, h/2};
Point(8) = {2.15, 0.3, 0, h/2};
Point(9) = {0, 1, 0, h};
Point(10) = {1, 1, 0, h/2};
Point(11) = {1.85, 1, 0, h};
Point(12) = {2.15, 1, 0, h};
Point(13) = {3, 1, 0, h/2};
Point(14) = {4, 1, 0, h};
Point(15) = {2, 0.3, 0, h/2};
Point(16) = {2, 1, 0, h/2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 7};

Line(5) = {8, 4};
Line(6) = {4, 5};
Line(7) = {5, 6};
Line(8) = {6, 14};
Line(9) = {14, 13};
Line(10) = {13, 12};

Line(12) = {11, 10};
Line(13) = {10, 9};
Line(14) = {9, 1};
Line(15) = {10, 2};
Line(16) = {2, 7};
Line(17) = {7, 10};
Line(18) = {7, 11};
Line(19) = {12, 8};
Line(20) = {8, 5};
Line(21) = {5, 13};
Line(22) = {13, 8};
Line(41) = {16, 15};
Line(42) = {7, 15};
Line(43) = {15, 8};
Line(44) = {12, 16};
Line(45) = {16, 11};

Line Loop(24) = {14, 1, -15, 13};
Plane Surface(24) = {24};
Line Loop(26) = {2, 3, -16};
Plane Surface(26) = {26};
Line Loop(28) = {16, 17, 15};
Plane Surface(28) = {28};
Line Loop(30) = {17, -12, -18};
Plane Surface(30) = {30};
Line Loop(34) = {22, -19, -10};
Plane Surface(34) = {34};
Line Loop(36) = {20, 21, 22};
Plane Surface(36) = {36};
Line Loop(38) = {5, 6, -20};
Plane Surface(38) = {38};
Line Loop(40) = {7, 8, 9, -21};
Plane Surface(40) = {40};

Line Loop(46) = {45, -18, 42, -41};
Plane Surface(47) = {46};
Line Loop(48) = {43, -19, 44, 41};
Plane Surface(49) = {48};


Physical Line("inlet") = {14};
Physical Line("gamma0") = {1, 7, 9, 13};
Physical Line("gamma26") = {2, 3, 5, 6};
Physical Line("gamma48") = {10, 12};
Physical Line("gamma5") = {42, 43, 44, 45};
Physical Line("inter23") = {16};
Physical Line("inter34") = {17};
Physical Line("inter45") = {18};
Physical Line("inter58") = {19};
Physical Line("inter78") = {22};
Physical Line("inter67") = {20};

Physical Surface("omega0") = {24, 40};
Physical Surface("omega26") = {26,38};
Physical Surface("omega3") = {28};
Physical Surface("omega48") = {30,34};
Physical Surface("omega5") = {47,49};
Physical Surface("omega7") = {36};


Physical Line("midflux") = {41};


