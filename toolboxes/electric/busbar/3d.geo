h=0.1;

Lx=1;
Ly=2;
Lz=0.5;

Point(1) = {0, -Ly/2.,-Lz/2.,h};
Point(2) = {0,  Ly/2.,-Lz/2.,h};
Point(3) = {Lx, Ly/2.,-Lz/2.,h};
Point(4) = {Lx,-Ly/2.,-Lz/2.,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
out[] = Extrude{0,0,Lz} {Surface{1};};


Physical Volume("omega") = {out[1]};
Physical Surface("top") = {out[0]};
Physical Surface("bottom") = {1};
Physical Surface("Rside") = {out[2]};
Physical Surface("V0") = {out[3]};
Physical Surface("Lside") = {out[4]};
Physical Surface("V1") = {out[5]};

