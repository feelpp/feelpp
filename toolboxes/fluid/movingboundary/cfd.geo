h=0.05;
Lenght=5;
Width=1;
Radius=2;
Thickness=0.1;
Point(1) = {0., 0., 0., h};
Point(2) = {Lenght, 0., 0., h};
Point(3) = {Lenght, Width, 0., h};
Point(4) = {0., Width, 0., h};

Point(9) = {Lenght, 0., 0., h};
Point(10) = {Lenght, Width, 0., h};

CenterX = Lenght+ Sqrt( -(Width/2)^2 + Radius^2) ;
Point(5) = {CenterX, Width/2., 0., h};
Point(6) = {CenterX, Width/2.+Radius, 0., h};
Point(7) = {CenterX, Width/2.-Radius, 0., h};
Point(8) = {CenterX+Radius, Width/2., 0., h};

Angle1Top=Pi/6.;
Angle2Top=Angle1Top+Pi/12.;
Angle1Bottom=-Pi/6.;
Angle2Bottom=Angle1Bottom-Pi/12.;

Point(11) = {CenterX+Radius*Cos(Angle1Top), Width/2.+Radius*Sin(Angle1Top), 0., h};
Point(12) = {CenterX+Radius*Cos(Angle2Top), Width/2.+Radius*Sin(Angle2Top), 0., h};
Point(13) = {CenterX+Radius*Cos(Angle1Top)+3, Width/2.+Radius*Sin(Angle1Top)+3, 0., h};
Point(14) = {CenterX+Radius*Cos(Angle2Top)+3, Width/2.+Radius*Sin(Angle2Top)+3, 0., h};

Point(15) = {CenterX+Radius*Cos(Angle1Bottom), Width/2.+Radius*Sin(Angle1Bottom), 0., h};
Point(16) = {CenterX+Radius*Cos(Angle2Bottom), Width/2.+Radius*Sin(Angle2Bottom), 0., h};
Point(17) = {CenterX+Radius*Cos(Angle1Bottom)+3, Width/2.+Radius*Sin(Angle1Bottom)-3, 0., h};
Point(18) = {CenterX+Radius*Cos(Angle2Bottom)+3, Width/2.+Radius*Sin(Angle2Bottom)-3, 0., h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {2, 5, 7};
Circle(6) = {7, 5, 16};
Circle(7) = {15, 5, 8};
Circle(8) = {8, 5, 11};
Circle(9) = {12, 5, 6};
Circle(10) = {6, 5, 3};
Line(11) = {12, 14};
Line(12) = {14, 13};
Line(13) = {13, 11};
Line(14) = {15, 17};
Line(15) = {17, 18};
Line(16) = {18, 16};
Circle(17) = {16, 5, 15};
Circle(18) = {11, 5, 12};
Line Loop(19) = {4, 1, 2, 3};
Plane Surface(20) = {19};
Line Loop(21) = {10, -2, 5, 6, 17, 7, 8, 18, 9};
Plane Surface(22) = {21};
Line Loop(23) = {18, 11, 12, 13};
Plane Surface(24) = {23};
Line Loop(25) = {17, 14, 15, 16};
Plane Surface(26) = {25};


Physical Line("wall_left") = {4};
Physical Line("outlet1") = {12};
Physical Line("outlet2") = {15};
Physical Line("wall_bottom") = {1};
Physical Line("wall_top") = {3};
Physical Line("wall_others") = {5,6,7,8,9,10,11,13,14,16};
Physical Surface("Fluid") = {20,22,24,26};
