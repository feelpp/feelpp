// A simple geometry for the lamina cribrosa and the central retinal
// vein passing through it

// The cgs unit system is adopted

// Outer radius of the lamina. Ref: Sigal 2005
LC_outer_radius = 0.095;

// Laminar thickness. Ref: Sigal 2005
LC_thickness    = 0.03;

// CRV radius. Ref: Jonas 1991
CRV_radius      = 0.008;

// Angle defining the portion of the LC lateral surface where we impose
// Neumann boundary conditions
alpha = Pi/3;

// Outer circle

Point(1) = {LC_outer_radius, 0, -LC_thickness/2};
Point(2) = {0, 0, -LC_thickness/2};
Point(3) = {LC_outer_radius*Cos(alpha), LC_outer_radius*Sin(alpha), -LC_thickness/2};
Point(4) = {LC_outer_radius*Cos(Pi+alpha/2), LC_outer_radius*Sin(Pi+alpha/2), -LC_thickness/2};

Circle(1) = {1, 2, 3};
Circle(2) = {3, 2, 4};
Circle(3) = {4, 2, 1};

// Inner hole

Point(5) = {CRV_radius, 0, -LC_thickness/2};
Point(6) = {0, CRV_radius, -LC_thickness/2};
Point(7) = {-CRV_radius, 0, -LC_thickness/2};
Point(8) = {0, -CRV_radius, -LC_thickness/2};

Circle(5) = {5, 2, 6};
Circle(6) = {6, 2, 7};
Circle(7) = {7, 2, 8};
Circle(8) = {8, 2, 5};

Extrude{0, 0, LC_thickness} { Line{1, 2, 3}; Line{-5,-8,-7,-6}; }

// Create all the remaining surfaces
Line Loop(1) = {1, 2, 3};
Line Loop(2) = {-5, -8, -7, -6};
Plane Surface(1000) = {1, 2};

Line Loop(3) = {9, 13, 17};
Line Loop(4) = {21, 25, 29, 33};
Plane Surface(1001) = {3, 4};


Surface Loop(2000) = {12, 16, 20, 1000, 1001, 24, 28, 32, 36};
Volume(3000) = {2000};

// Neumann and Dirichlet BC

// Fraction of the outer border of the LC covered by the Circle of Zin and Haller
Physical Surface("outer_ZHcircle") = {16, 20};

// Fraction of the outer border of the LC without the circle of Zin-Haller
Physical Surface("outer_noZHcircle") = {12};

// Inner border facing the Central Retinal Vein
Physical Surface("CRV") = {24, 28, 32, 36};

Physical Surface("Top") = {1001};
Physical Surface("Bottom") = {1000};

Physical Volume(999) = {3000};

