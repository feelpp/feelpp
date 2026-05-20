RCircle = 1.;
RDom = 10;

h = 3;
lcCircle = h/10;
lcDom = h;

// The construction of the center circle
Center = 0;
Centery = 0;

Point(9) = {Center,Centery,0,lcCircle};
Point(10) = {Center+RCircle,Centery,0,lcCircle};
Point(11) = {Center-RCircle,Centery,0,lcCircle};
Circle(7) = {10,9,11};
Circle(8) = {11,9,10};

// The construction of the left circle
leftCenter = Center-10*RCircle;

Point(1) = {leftCenter,Centery,0,lcCircle};
Point(2) = {leftCenter+RCircle,Centery,0,lcCircle};
Point(3) = {leftCenter,RCircle+Centery,0,lcCircle};
Point(4) = {leftCenter,-RCircle+Centery,0,lcCircle};
Circle(1) = {2,1,3};
Circle(2) = {4,1,2};
Circle(3) = {3,1,4};

RightCenter = Center+10*RCircle;

// The construction of the right circle
Point(5) = {RightCenter,Centery,0,lcCircle};
Point(6) = {RightCenter-RCircle,Centery,0,lcCircle};
Point(7) = {RightCenter,RCircle+Centery,0,lcCircle};
Point(8) = {RightCenter,-RCircle+Centery,0,lcCircle};
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

// Fluid domain
HalfSide = 50;
HalfHeight = 20;

Point(16) = {0,HalfHeight,0,lcDom};
Point(17) = {0,-HalfHeight,0,lcDom};
Point(18) = {0,HalfSide,0,lcDom};
Point(19) = {0,2*HalfSide+HalfHeight,0,lcDom};
Point(20) = {0,2*HalfSide-HalfHeight,0,lcDom};
Circle(14) = {16, 18, 20};
Circle(15) = {20, 18, 16};
Circle(16) = {17, 18, 19};
Circle(17) = {19, 18, 17};
Curve Loop(12) = {17, 16};
Curve Loop(13) = {15, 14};
Curve Loop(14) = {12, -12, 13,-13};
Plane Surface(12) = {7, 9, 11, 12, 13,14};


Physical Curve("CircleLeft") = {1, 3, 2};
Physical Curve("CircleCenter") = {7, 8};
Physical Curve("CircleRight") = {4, 5, 6};
Physical Curve("BoxWalls") = {17, 16, 15, 14};

Physical Surface("CirLeft") = {8};
Physical Surface("CirCenter") = {11};
Physical Surface("CirRight") = {10};
Physical Surface("Fluid") = {12};
