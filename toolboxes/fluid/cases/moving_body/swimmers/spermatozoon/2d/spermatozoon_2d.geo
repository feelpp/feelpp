// Channel with 7.5 micrometer height
Point(1) = {-20, 0, 0, 1.0};
Point(2) = {105, 0, 0, 1.0};
Point(3) = {105, 15, 0, 1.0};
Point(4) = {-20, 15, 0, 1.0};
// Channel with 15 micrometer height
//Point(1) = {-10, -7.5, 0, 1.0};
//Point(2) = {85, -7.5, 0, 1.0};
//Point(3) = {85, 22.5, 0, 1.0};
//Point(4) = {-10, 22.5, 0, 1.0};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Point(5) = {18, 7.5, -0, 0.5};
Point(6) = {16, 7.5, -0, 0.5};
Point(7) = {20, 7.5, -0, 0.5};
Point(8) = {18, 6.5, -0, 0.5};
Point(9) = {18, 8.5, -0, 0.5};
Point(11) = {20, 7.6, -0, 0.5};
Point(12) = {20, 7.4, -0, 0.5};

Ellipse(5) = {6, 5, 11, 9};
Ellipse(6) = {9, 5, 6, 11};
Ellipse(7) = {12, 5, 6, 8};
Ellipse(8) = {8, 5, 12, 6};

Point(13) = {70, 7.6, -0, 0.05};
Point(14) = {70, 7.4, -0, 0.05};
Line(9) = {11, 13};
Line(10) = {13, 14};
Line(11) = {14, 12};
Line Loop(12) = {4, 1, 2, 3};
Line Loop(13) = {5, 6, 9, 10, 11, 7, 8};
Plane Surface(14) = {12, 13};
Plane Surface(15) = {13};
Physical Line("Inlet") = {1};
Physical Line("Walls") = {4, 2};
Physical Line("Outlet") = {3};
Physical Line("Tail") = {9, 11, 10};
Physical Line("Head") = {7, 8, 5, 6};
Physical Surface("Fluid") = {14};
Physical Surface("Swimmer") = {15};
