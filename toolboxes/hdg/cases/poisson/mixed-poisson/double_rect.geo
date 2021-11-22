cl__1 = 1;

// h=0.25
// Mesh.CharacteristicLengthFactor=h

Point(1) = {0, 0, 0, 1};
Point(2) = {0, 1, 0, 1};
Point(4) = {2, 1, 0, 1};
Point(5) = {2, 0, 0, 1};
Line(1) = {2, 1};
Line(2) = {1, 5};
Line(3) = {5, 4};
Line(4) = {4, 2};
Line Loop(7) = {1, 2, 3, 4};
Plane Surface(8) = {7};
Extrude {0, 0, 1} {
  Surface{8};
}

// Second rectangle
Extrude {0, 0, 1} {
  Surface{30};
}

// Physical labels
Physical Surface("dirichlet1") = {21, 8, 25, 29, 17};
Physical Surface("interface") = {30};
Physical Surface("dirichlet2") = {43, 47, 39, 51, 52};
Surface Loop(53) = {43, 47, 51, 39, 52, 30};
Volume(54) = {53};
Physical Volume("omega1") = {1};
Surface Loop(55) = {21, 8, 17, 29, 25, 30};
Volume(56) = {55};
Physical Volume("omega2") = {2};
Physical Volume("omega") = {2, 1};

