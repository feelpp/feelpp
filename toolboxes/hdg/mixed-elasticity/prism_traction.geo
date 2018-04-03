h = 0.2;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {0, 1, 0, h};
Point(4) = {1, 1, 0, h};

Line(1) = {3, 1};
Line(2) = {1, 2};
Line(3) = {2, 4};
Line(4) = {4, 3};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Extrude {0, 0, 5} {
  Surface{6};
}


Physical Surface("bottom") = {6};
Physical Surface("top") = {28};
Physical Surface("lateral") = {15, 27, 23, 19};
Physical Volume("omega") = {1};
