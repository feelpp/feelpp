h = 0.1;

Point(1) = {-1,-1,-1,h};
Point(2) = {1,-1,-1,h};
Point(3) = {1,1,-1,h};
Point(4) = {-1,1,-1,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

out[] = Extrude{0,0,2} {Surface{1};};

Physical Volume("omega") = {out[1]};
Physical Surface("bottom") = {1};
Physical Surface("top") = {out[0]};
Physical Surface("side1") = {out[2]};
Physical Surface("side2") = {out[3]};
Physical Surface("side3") = {out[4]};
Physical Surface("side4") = {out[5]};
