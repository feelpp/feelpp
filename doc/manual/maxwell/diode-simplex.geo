h=0.25;
Point(1) = {0, 0, 0, h};
Point(2) = {0, 1, 0, h};
Point(3) = {0, 2, 0, h};
Point(4) = {2, 0, 0, h};
Point(5) = {2, 1, 0, h};
Point(6) = {2, 2, 0, h};
Point(7) = {3, 0, 0, h};
Point(8) = {4, 0, 0, h};
Point(9) = {4, 2, 0, h};
Point(10) = {2.707106781, 0.707106781, 0, h};

Line(1) = {2, 3};
Line(2) = {3, 9};
Line(3) = {9, 8};
Line(4) = {8, 7};
Circle(5) = {7, 4, 5};
Line(6) = {5, 2};

Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

Physical Line("Dirichlet") = {1}
Physical Line("Metal") = {2,3,4,5,6}

