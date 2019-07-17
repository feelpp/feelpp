h = 0.1;

Point(1) = {0, 0, 0, h};
Point(2) = {0, 0, 5, h};
Point(3) = {-1, 0, 5, h};
Point(4) = {-1, 0, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(5) = {5};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{5};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{22};
}


Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{39};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{56};
}
Physical Line("centerline") = {1};
Physical Surface("inlet") = {21, 38, 55, 72};
Physical Surface("outlet") = {14, 31, 48, 65};
Physical Surface("wall") = {18, 35, 52, 69};
Physical Volume("wall") = {1, 2, 3, 4};
