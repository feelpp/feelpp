h=0.03;
Point(1) = {0., 0., 0., h};
Point(2) = {0, 1., 0., h};
Point(3) = {1, 1, 0., h};
Point(4) = {1., 0, 0., h};

Point(5) = {0.5, 0, 0., h};
Point(6) = {0.5, 1., 0., h};
Point(7) = {0., 0.5, 0., h};
Point(8) = {1., 0.5, 0., h};
Point(9) = {0.5, 0.5, 0., h};

Line(1) = {2, 7};
Line(2) = {7, 1};
Line(3) = {1, 5};
Line(4) = {5, 4};
Line(5) = {4, 8};
Line(6) = {8, 3};
Line(7) = {3, 6};
Line(8) = {6, 2};
Line(9) = {6, 9};
Line(10) = {9, 5};
Line(11) = {7, 9};
Line(12) = {9, 8};
Curve Loop(1) = {8, 1, 11, -9};
Plane Surface(1) = {1};
Curve Loop(2) = {9, 12, 6, 7};
Plane Surface(2) = {2};
Curve Loop(3) = {12, -5, -4, -10};
Plane Surface(3) = {3};
Curve Loop(4) = {10, -3, -2, 11};
Plane Surface(4) = {4};

Physical Line("insulated") = {3,4,7,8};
Physical Line("hot_wall") = {1,2};
Physical Line("cold_wall") = {5,6};
Physical Line("horizontal_centerline") = {11,12};
Physical Line("vertical_centerline") = {9,10};

Physical Surface("Omega") = {1,2,3,4};
