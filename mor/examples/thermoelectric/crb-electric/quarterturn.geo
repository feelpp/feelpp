// Define Main params
Unit = 1.e-3;
lc = 1*Unit;
lc_ext = 3*lc;
lc_inf = 1*lc;

h=0.2;
r1=1;
r2=2;
L=2*r2;



//Mesh.ElementOrder = 1;
Point(1) = {0, 0, -L, h};
Point(2) = {r1, 0, -L, h};
Point(3) = {r2, 0, -L, h};
Point(4) = {0, r1, -L, h};
Point(5) = {0, r2, -L, h};
Circle(1) = {2, 1, 4};
Circle(2) = {3, 1, 5};
Line(3) = {4, 5};
Line(4) = {2, 3};
Line Loop(5) = {3, -2, -4, 1};
Plane Surface(1) = {5};

out[] = Extrude {0,0,L} {Surface{1};};

Physical Volume("omega") = {out[1]};
Physical Surface("top") = {out[0]};
Physical Surface("bottom") = {1};
Physical Surface("Rint") = {out[5]};
Physical Surface("Rext") = {out[3]};
Physical Surface("in") = {out[2]};
Physical Surface("out") = {out[4]};
