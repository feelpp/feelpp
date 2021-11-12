// Gmsh project created on Tue Oct  5 09:43:44 2021
SetFactory("OpenCASCADE");

h = 0.0075;

r1 = 1.9;
r2 = 2.1;
R = r2-r1;
angle = Pi/6;
startangle = Pi/2;
endangle = startangle + angle;

Point(1) = {0, 0, 0, h};
Point(2) = {r1*Cos(startangle), r1*Sin(startangle), 0, h};
Point(3) = {r2*Cos(startangle), r2*Sin(startangle), 0, h};
Point(4) = {r2*Cos(endangle), r2*Sin(endangle), 0, h};
Point(5) = {r1*Cos(endangle), r1*Sin(endangle), 0, h};

Line(1) = {2, 3};
Line(2) = {5, 4};
Circle(3) = {2, 1, 5};
Circle(4) = {3, 1, 4};
Line Loop(1) = {3, 2, -4, -1};
Surface(1) = {1};

Physical Line("wall") = {4, 3};
Physical Line("outlet") = {1};
Physical Line("inlet") = {2};
Physical Surface("Omega") = {1};