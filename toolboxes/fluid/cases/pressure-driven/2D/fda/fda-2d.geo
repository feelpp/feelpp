// Gmsh project created on Fri Oct  1 15:05:58 2021
SetFactory("OpenCASCADE");
h = 0.001;
di = 0.012;
t = 0.004;
cc = 0.022685;
Th = 0.04;
Li = 0.1;
Lo = 0.2;

Point(1) = {0, 0, 0, h};
Point(2) = {0, di, 0, h};
Point(3) = {Li, 0, 0, h};
Point(4) = {Li, di, 0, h};
Point(5) = {Li+cc, 0.5*(di-t), 0, h};
Point(6) = {Li+cc, 0.5*(di+t), 0, h};
Point(7) = {Li+cc+Th, 0.5*(di-t), 0, h};
Point(8) = {Li+cc+Th, 0.5*(di+t), 0, h};
Point(9) = {Li+cc+Th, 0, 0, h};
Point(10) = {Li+cc+Th, di, 0, h};
Point(11) = {Li+cc+Th+Lo, 0, 0, h};
Point(12) = {Li+cc+Th+Lo, di, 0, h};

Line(1) = {2, 1};
Line(2) = {12, 11};
Line(3) = {1, 3};
Line(4) = {3, 5};
Line(5) = {5, 7};
Line(6) = {7, 9};
Line(7) = {9, 11};
Line(8) = {12, 10};
Line(9) = {10, 8};
Line(10) = {8, 6};
Line(11) = {6, 4};
Line(12) = {4, 2};

Line Loop(1) = {12, 1, 3, 4, 5, 6, 7, -2, 8, 9, 10, 11};

Plane Surface(1) = {1};

Physical Line("inlet") = {1};
Physical Line("outlet") = {2};
Physical Line("wall") = {12, 11, 10, 9, 8, 3, 4, 5, 6, 7};

Physical Surface("Omega") = {1};
