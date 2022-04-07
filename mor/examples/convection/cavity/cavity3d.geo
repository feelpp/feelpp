h = 0.025;
Point(1) = {0,0,0.0,h};
Point(2) = {1,0,0.0,h};
Point(3) = {1,1,0.0,h};
Point(4) = {0,1,0.0,h};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(7) = {1,2,3,4};
Plane Surface(8) = {7};

Extrude {0, 0, 1} {Surface{8};}
Physical Surface("wallInsulated") = {30,17,8,25};
Physical Surface("wallRight") = {21};
Physical Surface("wallLeft") = {29};
Physical Volume("domain") = {1};

