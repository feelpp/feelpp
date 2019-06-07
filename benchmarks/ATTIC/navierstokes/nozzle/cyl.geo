h=0.05;
Point(1) = {0, 0, 0, h};
Point(2) = {0, 1, 0, h};
Point(3) = {5, 1, 0, h};
Point(4) = {5, 0, 0, h};
Point(5) = {1.15, 0.6, 0, h};
Point(6) = {1.2, 0.6, 0, h};
Point(7) = {1.25, 0.6, 0, h};

Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};


Circle(5) = {7, 6, 5};
Circle(6) = {5, 6, 7};
Line Loop(7) = {3, 4, 1, 2};
Line Loop(8) = {6, 5};
Plane Surface(9) = {7, 8};
Physical Line("inlet") = {4};
Physical Line("wall") = {1, 3, 5, 6};
Physical Line("outlet") = {2};
Physical Surface("fluid") = {9};
