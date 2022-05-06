h=0.1;

Lx=1;
Ly=2;

Point(1) = {0, -Ly/2.,0,h};
Point(2) = {0,  Ly/2.,0,h};
Point(3) = {Lx, Ly/2.,0,h};
Point(4) = {Lx,-Ly/2.,0,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Physical Surface("omega") = {1};
Physical Line("V0") = {1};
Physical Line("Lside") = {2};
Physical Line("V1") = {3};
Physical Line("Rside") = {4};

