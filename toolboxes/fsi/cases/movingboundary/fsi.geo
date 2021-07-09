h=0.05;
Width=1;
Radius=2;
Thickness=0.1;
Lenght=5+Thickness;
Point(1) = {0., 0., 0., h};
Point(2) = {Lenght-Thickness, 0., 0., h};
Point(3) = {Lenght-Thickness, Width, 0., h};
Point(4) = {0., Width, 0., h};

Point(9) = {Lenght, 0., 0., h};
Point(10) = {Lenght, Width, 0., h};

CenterX = Lenght+ Sqrt( -(Width/2)^2 + Radius^2) ;
Point(5) = {CenterX, Width/2., 0., h};
Point(6) = {CenterX, Width/2.+Radius, 0., h};
Point(7) = {CenterX, Width/2.-Radius, 0., h};
Point(8) = {CenterX+Radius, Width/2., 0., h};

Point(11) = {CenterX, Width/2.+Radius+Thickness, 0., h};
Point(12) = {CenterX, Width/2.-Radius-Thickness, 0., h};
Point(13) = {CenterX+Radius+Thickness, Width/2., 0., h};

Line(1) = {1, 2};
Line(2) = {2, 9};
Line(3) = {9, 10};
Line(4) = {10, 3};
Line(5) = {3, 4};
Line(6) = {4, 1};
Circle(7) = {9, 5, 7};
Circle(8) = {7, 5, 8};
Circle(9) = {8, 5, 6};
Circle(10) = {6, 5, 10};
Line Loop(1) = {6, 1, 2, 3, 4, 5};
Plane Surface(1) = {1};
Line Loop(2) = {10, -3, 7, 8, 9};
Plane Surface(2) = {2};

Circle(11) = {3, 5, 11};
Circle(12) = {11, 5, 13};
Circle(13) = {13, 5, 12};
Circle(14) = {12, 5, 2};
Line Loop(3) = {11, 12, 13, 14, 2, 7, 8, 9, 10, 4};
Plane Surface(3) = {3};

Physical Line("wall_left") = {6};
Physical Line("wall_bottom") = {1};
Physical Line("wall_top") = {5};
Physical Line("wall_balloon_internal") = {7,8,9,10};
Physical Line("wall_balloon_external") = {11,12,13,14};
Physical Line("wall_balloon_fixed") = {2,4};
Physical Surface("Fluid") = {1,2};
Physical Surface("Solid") = {3};
