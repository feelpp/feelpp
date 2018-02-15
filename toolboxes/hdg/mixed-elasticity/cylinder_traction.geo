h = 0.2;

Point(1) = {0, 0, 0, h};
Point(7) = {1, 0, 0, h};
Point(8) = {-1, 0, 0, h};
Point(9) = {0, 1, 0, h};
Point(10) = {0, -1, 0, h};
Circle(5) = {8, 1, 9};
Circle(6) = {9, 1, 7};
Circle(7) = {7, 1, 10};
Circle(8) = {10, 1, 8};

Line Loop(10) = {5, 6, 7, 8};
Plane Surface(12) = {10};


Plane Surface(13) = {10};
Extrude {0, 0, 5} {
  Surface{12};
}

Physical Volume("omega") = {1};
Physical Surface("lateral") = {34, 30, 26, 22};
Physical Surface("top") = {35};
Physical Surface("bottom") = {12};




