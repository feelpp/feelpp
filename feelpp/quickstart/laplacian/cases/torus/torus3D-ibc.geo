h = 0.08;
ri = 1.0;
re = 2.0;
H = 1.0;
delta = Pi/2.0;

Mesh.ElementOrder = 1;

Point(1) = {0, 0, 0, h};
Point(2) = {ri, 0, 0, h};
Point(3) = {re, 0, 0, h};
Point(4) = {re, 0, H, h};
Point(5) = {ri, 0, H, h};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Extrude {{0, 0, 1}, {0, 0, 0}, delta} {
  Surface{1};
  Layers{15};
}

Physical Volume("Omega") = {1};
Physical Surface("Dirichlet") = {1};
Physical Surface("Ibc") = {26};
Physical Surface("Neumann") = {13, 17, 21, 25};
