// Mesh size
h=0.025;

// Boundaries
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 2, 0, h};
Point(4) = {0, 2, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Physical Line("Bottom") = {1};
Physical Line("Right") = {2};
Physical Line("Top") = {3};
Physical Line("Left") = {4};

// Fluid surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface("OmegaFluid") = {1};

