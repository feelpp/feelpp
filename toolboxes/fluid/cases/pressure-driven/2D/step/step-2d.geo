SetFactory("OpenCASCADE");

h=0.1;

hin=1;
hout=2;
lstep=1;
lfloor=3;

Point(1) = {-lstep, hout-hin, 0, h};
Point(2) = {0, hout-hin, 0, h*0.1};
Point(3) = {0, 0, 0, h*0.5};
Point(4) = {lfloor, 0, 0, h};
Point(5) = {lfloor, hout, 0, h};
Point(6) = {-lstep, hout, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

Physical Line("wall") = {1, 2, 3, 5};
Physical Line("inlet") = {6};
Physical Line("outlet") = {4};
Physical Surface("Omega") = {1};