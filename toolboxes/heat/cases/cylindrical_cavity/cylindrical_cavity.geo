// Simple closed cylindrical cavity benchmark.
//
// The geometry is intentionally built with the built-in kernel and extrusion
// return values. This avoids depending on OpenCASCADE entity numbers after
// BooleanFragments, which changed across Gmsh versions.

h = 0.1;

R = 0.5;
Rext = 0.6;
L = 1.0;
Z = 0.1;

// Bottom conducting disk seed surface: z = -Z.
Point(1) = {0, 0, -Z, h};
Point(2) = {R, 0, -Z, h};
Point(3) = {0, R, -Z, h};
Point(4) = {-R, 0, -Z, h};
Point(5) = {0, -R, -Z, h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Top conducting disk seed surface: z = L + Z.
Point(101) = {0, 0, L + Z, h};
Point(102) = {R, 0, L + Z, h};
Point(103) = {0, R, L + Z, h};
Point(104) = {-R, 0, L + Z, h};
Point(105) = {0, -R, L + Z, h};
Circle(101) = {102, 101, 103};
Circle(102) = {103, 101, 104};
Circle(103) = {104, 101, 105};
Circle(104) = {105, 101, 102};
Curve Loop(101) = {101, 102, 103, 104};
Plane Surface(101) = {101};

// Lateral conducting shell seed annulus: R <= r <= Rext, z = 0.
Point(201) = {0, 0, 0, h};
Point(202) = {Rext, 0, 0, h};
Point(203) = {0, Rext, 0, h};
Point(204) = {-Rext, 0, 0, h};
Point(205) = {0, -Rext, 0, h};
Point(206) = {R, 0, 0, h};
Point(207) = {0, R, 0, h};
Point(208) = {-R, 0, 0, h};
Point(209) = {0, -R, 0, h};

Circle(201) = {202, 201, 203};
Circle(202) = {203, 201, 204};
Circle(203) = {204, 201, 205};
Circle(204) = {205, 201, 202};
Circle(205) = {206, 201, 207};
Circle(206) = {207, 201, 208};
Circle(207) = {208, 201, 209};
Circle(208) = {209, 201, 206};

Curve Loop(201) = {201, 202, 203, 204};
Curve Loop(202) = {-208, -207, -206, -205};
Plane Surface(201) = {201, 202};

bottom[] = Extrude {0, 0, Z} {
  Surface{1};
};

top[] = Extrude {0, 0, -Z} {
  Surface{101};
};

side[] = Extrude {0, 0, L} {
  Surface{201};
};

Physical Volume("MaterialBottom", 69) = {bottom[1]};
Physical Volume("MaterialTop", 407) = {top[1]};
Physical Volume("LateralVolume", 606) = {side[1]};

Physical Surface("CavityBottom", 71) = {bottom[0]};
Physical Surface("CavityTop", 409) = {top[0]};
Physical Surface("CavitySide", 410) = {side[6], side[7], side[8], side[9]};

Physical Surface("ExternalBoundaryBottom", 72) = {1};
Physical Surface("ExternalBoundaryTop", 607) = {101};
Physical Surface("ExternalBoundaryLateral", 608) = {side[2], side[3], side[4], side[5]};

Characteristic Length { PointsOf{ Surface{:}; } } = h;
