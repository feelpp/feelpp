h=0.1;

Point(1) = {0,0,0,h};
Point(2) = {1,0,0,h};
Point(3) = {0,1,0,h};
Point(4) = {0,2,0,h};
Point(5) = {2,0,0,h};

Circle(1) = {2,1,3};
Line(2) = {3,4};
Circle(3) = {4,1,5};
Line(4) = {5,2};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Physical Surface("omega") = {1};
Physical Line("Rint") = {1};
Physical Line("V0") = {2};
Physical Line("Rext") = {3};
Physical Line("V1") = {4};

