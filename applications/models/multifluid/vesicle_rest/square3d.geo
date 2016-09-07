h = 0.1;
xmin = -4;
xmax = 4;
ymin = -4;
ymax = 4;
zmin = -4;
zmax = 4;
Point(1) = {xmin,ymin,zmin,h};
Point(2) = {xmax,ymin,zmin,h};
Point(3) = {xmax,ymax,zmin,h};
Point(4) = {xmin,ymax,zmin,h};
Point(5) = {xmin,ymin,zmax,h};
Point(6) = {xmax,ymin,zmax,h};
Point(7) = {xmax,ymax,zmax,h};
Point(8) = {xmin,ymax,zmax,h};
Line(1) = {4,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {8,5};
Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {7,8};
Line(9) = {4,8};
Line(10) = {1,5};
Line(11) = {2,6};
Line(12) = {3,7};

Line Loop(101) = {2,1,4,3};
Line Loop(102) = {7,11,-3,-12};
Line Loop(103) = {5,6,7,8};
Line Loop(104) = {1,10,-5,-9};
Line Loop(105) = {10,6,-11,-2};
Line Loop(106) = {4,9,-8,-12};

Plane Surface(201) = {101};
Plane Surface(202) = {102};
Plane Surface(203) = {103};
Plane Surface(204) = {104};
Plane Surface(205) = {105};
Plane Surface(206) = {106};

Surface Loop(301) = {201, 202, 203, 204, 205, 206};

Physical Surface("Left") = {204};
Physical Surface("Bottom") = {205};
Physical Surface("Right") = {202};
Physical Surface("Top") = {206};
Physical Surface("Front") = {203};
Physical Surface("Back") = {201};

Physical Volume("Omega") = {301};
