h = 0.1;

Point(1) = {0, 0, 0, h};
Point(2) = {0, 0, 1, h};
Point(3) = {0, 0, 2, h};
Point(4) = {0, 0, 5, h};
Point(5) = {-1, 0, 5, h};
Point(6) = {-1, 0, 2, h};
Point(7) = {-1, 0, 1, h};
Point(8) = {-1, 0, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 7};
Line(3) = {7, 8};
Line(4) = {8, 1};
Line(5) = {2, 3};
Line(6) = {3, 6};
Line(7) = {6, 7};
Line(8) = {3, 4};
Line(9) = {4, 5};
Line(10) = {5, 6};
//Line(11) = {1, 2, 3, 4};
Line Loop(11) = {1, 2, 3, 4};
Plane Surface(11) = {11};
Line Loop(12) = {5, 6, 7, -2};
Plane Surface(12) = {12};
Line Loop(13) = {8, 9, 10, -6};
Plane Surface(13) = {13};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{11, 12, 13};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{30, 47, 64};
}


Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{81, 98, 115};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{132, 149, 166};
}
Physical Line("centerline") = {1,5,8};
Physical Surface("inlet") = {29, 182, 131, 80};
Physical Surface("section1") = {175, 124, 73, 22};
Physical Surface("section2") = {141, 90, 39, 191};
Physical Surface("outlet") = {158, 204, 56, 107};
Physical Surface("wall") = {162, 208, 60, 111, 145, 195, 43, 94, 128, 179, 26, 77};
Physical Volume("fluid") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

