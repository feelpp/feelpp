// Gmsh project created on Tue Oct  5 09:43:44 2021
SetFactory("OpenCASCADE");

h = 0.1;

r1 = 1.9;
r2 = 2.1;
R = r2-r1;
angle = Pi/6;

Point(1) = {0, 0, 0, h};
Point(2) = {r1, 0, 0, h};
Point(3) = {r2, 0, 0, h};
Point(4) = {r2*Cos(angle), r2*Sin(angle), 0, h};
Point(5) = {r1*Cos(angle), r1*Sin(angle), 0, h};

Line(1) = {2, 3};
Line(2) = {5, 4};
Circle(3) = {2, 1, 5};
Circle(4) = {3, 1, 4};
Line Loop(1) = {3, 2, -4, -1};
Surface(1) = {1};

Extrude {0, 0, R} {
  Surface{1}; 
}

Physical Surface("wall") = {4, 2, 6, 1};
Physical Surface("outlet") = {3};
Physical Surface("inlet") = {5};
Physical Volume("Omega") = {1};