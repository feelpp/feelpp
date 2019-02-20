Point(1) = {0, 0, 0, h};
//+
Point(2) = {1, 0, 0, h};
//+
Point(3) = {1, 1, 0, h};
//+
Point(4) = {0, 1, 0, h};
//+
Point(5) = {0, 0.4, 0, h};
//+
Point(6) = {0, 0.6, 0, h};
//+
Point(7) = {1, 0.6, 0, h};
//+
Point(8) = {1, 0.4, 0, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 8};
//+
Line(3) = {8, 5};
//+
Line(4) = {5, 1};
//+
Line(5) = {6, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 8};
//+
Line(9) = {7, 3};
//+
Line(10) = {3, 4};
//+
Line(11) = {4, 6};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {3, -6, 7, 8};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {9, 10, 11, 7};
//+
Plane Surface(3) = {3};
//+
Physical Line("gamma1") = {4, 1, 2};
//+
Physical Line("gamma2") = {11, 10, 9};
//+
Physical Line("gammainter") = {6, 8};
//+
Physical Line("interface1") = {3};
//+
Physical Line("interface2") = {7};
//+
Physical Surface("omega1") = {1};
//+
Physical Surface("omega2") = {3};
//+
Physical Surface("plate") = {2};