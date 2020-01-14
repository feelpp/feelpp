h=0.2;
r1=1;
r2=2;
Mesh.ElementOrder = 1;
Point(1) = {0, 0, 0, h};
Point(2) = {r1, 0, 0, h};
Point(3) = {r2, 0, 0, h};
Point(4) = {0, r1, 0, h};
Point(5) = {0, r2, 0, h};
Circle(1) = {2, 1, 4};
Circle(2) = {3, 1, 5};
Line(3) = {4, 5};
Line(4) = {2, 3};
Line Loop(5) = {3, -2, -4, 1};
Plane Surface(1) = {5};

out[] = Extrude {0,0,1} {Surface{1};};

Physical Volume("omega") = {out[1]};
Physical Surface("bottom") = {out[0]};
Physical Surface("top") = {1};
Physical Surface("R") = {out[3],out[5]};
Physical Surface("in") = {out[2]};
Physical Surface("out") = {out[4]};
