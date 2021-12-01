// Gmsh project created on Tue Oct  5 09:43:44 2021
SetFactory("OpenCASCADE");
SetFactory("OpenCASCADE");

h = 0.0075;

r1 = 1.9;
r2 = 2.1;
R = r2-r1;
angle = Pi/6;
alpha1 = Pi/2;
alpha2 = alpha1 + angle;

Point(1) = {0, 0, 0, h};
Point(2) = {r1*Cos(alpha1), r1*Sin(alpha1), 0, h};
Point(3) = {r2*Cos(alpha1), r2*Sin(alpha1), 0, h};
Point(4) = {r2*Cos(alpha2), r2*Sin(alpha2), 0, h};
Point(5) = {r1*Cos(alpha2), r1*Sin(alpha2), 0, h};

Line(1) = {2, 3};
Line(2) = {5, 4};
Circle(3) = {2, 1, 5};
Circle(4) = {3, 1, 4};
Line Loop(1) = {3, 2, -4, -1};
Surface(1) = {1};

Characteristic Length{ PointsOf{ Surface{1}; } } = h;

Extrude {0, 0, R} {
  Surface{1}; 
}

Physical Surface("wall") = {5, 3, 1, 6};
Physical Surface("outlet") = {2};
Physical Surface("inlet") = {4};
Physical Volume("Omega") = {1};