RCircle = 0.125;
h = 0.01;
lcCircle = 0.01;
lcDom = h;

// The construction of the center circle
Centerx = 3;
Centery = 0.29;

Point(9) = {Centerx,Centery,0,lcCircle};
Point(10) = {Centerx+RCircle,Centery,0,lcCircle};
Point(11) = {Centerx-RCircle,Centery,0,lcCircle};
Circle(7) = {10,9,11};
Circle(8) = {11,9,10};

// The construction of the left circle
theta1 = 23*Pi/24;
Centerx2 = Centerx + 10*RCircle*Cos(theta1);
Centery2 = Centery + 10*RCircle*Sin(theta1);

Point(1) = {Centerx2,Centery2,0,lcCircle};
Point(2) = {Centerx2+RCircle,Centery2,0,lcCircle};
Point(3) = {Centerx2,RCircle+Centery2,0,lcCircle};
Point(4) = {Centerx2,-RCircle+Centery2,0,lcCircle};
Circle(1) = {2,1,3};
Circle(2) = {4,1,2};
Circle(3) = {3,1,4};

// The construction of the right circle

theta2 = -Pi/24;
Centerx3 = Centerx + 10*RCircle*Cos(theta2);
Centery3 = Centery + 10*RCircle*Sin(theta2);

Point(5) = {Centerx3,Centery3,0,lcCircle};
Point(6) = {Centerx3-RCircle,Centery3,0,lcCircle};
Point(7) = {Centerx3,RCircle+Centery3,0,lcCircle};
Point(8) = {Centerx3,-RCircle+Centery3,0,lcCircle};
Circle(4) = {7,5,6};
Circle(5) = {6,5,8};
Circle(6) = {8,5,7};

// Defining the arcs and surfaces of the three circles
Line Loop(7) = {1,2,3};
Plane Surface(8) = {7};
Line Loop(9) = {4,5,6};
Plane Surface(10) = {9};
Line Loop(11) = {7,8};
Plane Surface(11) = {11};

// Constructing the two arms linking the three spheres
Line(12) = {11, 2};
Line(13) = {10, 6};

// Rectangle vertices
Point(12) = {0,0,0,lcDom};
Point(13) = {6,0,0,lcDom};
Point(14) = {6,2.,0,lcDom};
Point(15) = {0,2.,0,lcDom};

// Rectangle lines
Line(16) = {12, 13};
Line(17) = {13, 14};
Line(18) = {14, 15};
Line(19) = {15, 12};
Line Loop(20) = {16,17,18,19};

// Defining the rectangle surface
Plane Surface(21) = {20, 7, 9, 11};


Physical Curve("CircleLeft") = {1, 3, 2};
Physical Curve("CircleCenter") = {7, 8};
Physical Curve("CircleRight") = {4, 5, 6};
Physical Curve("BoxWalls") = {16, 17, 18, 19};

Physical Surface("CirLeft") = {8};
Physical Surface("CirCenter") = {11};
Physical Surface("CirRight") = {10};
Physical Surface("Fluid") = {21};