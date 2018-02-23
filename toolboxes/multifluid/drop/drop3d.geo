// Mesh size
h=0.05;

// Cylinder boundaries
Circle_xc = 0.5;
Circle_yc = 0.5;
Circle_r = 0.5;
Cyl_height = 2;

Point(1) = {Circle_xc, Circle_yc, 0, h};
Point(2) = {Circle_xc, Circle_yc-Circle_r, 0, h};
Point(3) = {Circle_xc+Circle_r, Circle_yc, 0, h};
Point(4) = {Circle_xc, Circle_yc+Circle_r, 0, h};
Point(5) = {Circle_xc-Circle_r, Circle_yc, 0, h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

cyl[] = Extrude {0, 0, Cyl_height} { Surface{1}; };

Physical Surface("Top") = { cyl[0] };
Physical Surface("Bottom") = {1};
Physical Surface("Side") = { cyl[2], cyl[3], cyl[4], cyl[5] };

// volume fluid
Physical Volume("OmegaFluid") = { cyl[1] };
