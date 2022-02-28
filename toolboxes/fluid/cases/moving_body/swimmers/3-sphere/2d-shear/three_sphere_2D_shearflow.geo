RCircle = 1.;
RDom = 10;

h = 3;
lcCircle = h/10;
lcDom = h;

// The construction of the center circle
Center = 0;
Centery = 0;
Beta= 30/180*3.141592;

Point(9) = {Center,Centery,0,lcCircle};
Point(10) = {Center+RCircle,Centery,0,lcCircle};
Point(11) = {Center-RCircle,Centery,0,lcCircle};
Circle(7) = {10,9,11};
Circle(8) = {11,9,10};

// The construction of the left circle
leftCenter = Center-10*RCircle;
CenterX_left=-15*Cos(Beta);
CenterY_left=-15*Sin(Beta);

Point(1) = {CenterX_left, CenterY_left,0,lcCircle};
Point(2) = {CenterX_left +RCircle, CenterY_left,0,lcCircle};
Point(3) = {CenterX_left,RCircle+ CenterY_left,0,lcCircle};
Point(4) = {CenterX_left,-RCircle+ CenterY_left,0,lcCircle};
Circle(1) = {2,1,3};
Circle(2) = {4,1,2};
Circle(3) = {3,1,4};

CenterX_right=15*Cos(Beta);
CenterY_right=15*Sin(Beta);


// The construction of the right circle
Point(5) = {CenterX_right, CenterY_right,0,lcCircle};
Point(6) = {CenterX_right-RCircle, CenterY_right,0,lcCircle};
Point(7) = {CenterX_right,RCircle+ CenterY_right,0,lcCircle};
Point(8) = {CenterX_right,-RCircle+ CenterY_right,0,lcCircle};

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


// Rectangle parameters

HalfSide = 30;
HalfHeight = 20;
//HalfSide = 5*RightCenter;
//HalfHeight = 2*RightCenter;

// Rectangle vertices
Point(12) = {-HalfSide,-HalfHeight,0,lcDom};
Point(13) = {HalfSide,-HalfHeight,0,lcDom};
Point(14) = {HalfSide,HalfHeight,0,lcDom};
Point(15) = {-HalfSide,HalfHeight,0,lcDom};

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
Physical Curve("Walls") = { 17,19};
Physical Curve("Bottom") = {16};
Physical Curve("Top") = {18};

Physical Surface("CirLeft") = {8};
Physical Surface("CirCenter") = {11};
Physical Surface("CirRight") = {10};
Physical Surface("Fluid") = {21};

