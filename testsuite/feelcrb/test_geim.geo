h = 0.1;

Point(1) = {0,0,0,h};
Point(2) = {0,1,0,h};
Point(3) = {1,1,0,h};
Point(4) = {1,0,0,h};

Point(5) = {0.5,0.5,0,h};
Point(6) = {0.8,0.5,0,h};
Point(7) = {0.5,0.8,0,h};
Point(8) = {0.2,0.5,0,h};
Point(9) = {0.5,0.2,0,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};

Surface(1) = {1,2};
Surface(2) = {2};

Physical Surface("Omega1") = {1};
Physical Surface("Omega2") = {2};

Physical Line("left") = {1};
Physical Line("top") = {2};
Physical Line("right") = {3};
Physical Line("bottom") = {4};
