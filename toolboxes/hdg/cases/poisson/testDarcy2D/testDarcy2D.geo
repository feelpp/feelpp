h=0.1;

Point(1) = {0,0,0,h};
Point(2) = {0,2,0,h};
Point(3) = {2,2,0,h};
Point(4) = {2,0,0,h};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};
Plane Surface(1) = {5};
Physical Surface("omega") = {1};
Physical Line("left") = {1};
Physical Line("top") = {2};
Physical Line("right") = {3};
Physical Line("bottom") = {4};
