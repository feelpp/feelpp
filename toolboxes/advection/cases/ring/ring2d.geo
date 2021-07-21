h=0.5;
Point(1) = {0.1,0,0,h/2};
Point(2) = {1,0,0,h};
Point(3) = {0,1,0,h};
Point(4) = {0,0.1,0,h/2};
Point(5) = {0,0,0,h/2};
Line(1) = {1,2};
Circle(2) = {2,5,3};
Line(3) = {3,4};
Circle(4) = {4,5,1};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Physical Line("Bottom") = {1};
Physical Line("OuterCircle") = {2};
Physical Line("Left") = {3};
Physical Line("InnerCircle") = {4};
Physical Surface("Omega") = {6};
