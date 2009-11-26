
hfine=5e-4;
h=1e-3;
Point(1) = {0,0,0,hfine};
Point(2) = {0.002,0,0,hfine};
Point(3) = {0.002,0.02,0,hfine};
Point(4) = {0.01,0.02,0,h};
Point(5) = {0.01,0.138,0,h};
Point(6) = {0.002,0.138,0,hfine};
Point(7) = {0.002,0.188,0,h};
Point(8) = {0,0.188,0,hfine};
Point(9) = {0,0.0138,0,hfine};
Point(10) = {0,0.0138,0,hfine};
Point(11) = {0,0.138,0,hfine};
Point(12) = {0,0.02,0,hfine};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,11};
Line(9) = {11,12};
Line(10) = {12,1};
Line(11) = {3,12};
Line(12) = {6,11};
Line Loop(13) = {1,2,11,10};
Plane Surface(14) = {13};
Line Loop(15) = {3,4,5,12,9,-11};
Plane Surface(16) = {15};
Line Loop(17) = {8,-12,6,7};
Plane Surface(18) = {17};


Physical Line("inflow") = {1};
Physical Line("wall") = {2,3,4,5,6};
Physical Line("outflow") = {7};
Physical Line("symmetry") = {8,9,10};

Physical Surface("dom1") = {14};
Physical Surface("dom2") = {16};
Physical Surface("dom3") = {18};

