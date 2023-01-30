A = 0.2; // large axis
B = 0.1; // small axis
h=0.01;
//+
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 0.6, 0, h};
Point(4) = {0, 0.6, 0, h};

Point(5) = {0.5, 0.1, 0, h};
Point(6) = {0.5-A*0.5, 0.1, 0, h};
Point(7) = {0.5+A*0.5, 0.1, 0, h};
Point(8) = {0.5, 0.1+0.5*B, 0, h};
Point(9) = {0.5, 0.1-0.5*B, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {4, 1, 2, 3};

Ellipse(5) = {7, 5, 8};
Ellipse(6) = {8, 5, 6};
Ellipse(7) = {6, 5, 9};
Ellipse(8) = {9, 5, 7};

Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {2};
Plane Surface(2) = {1, 2};

Physical Curve("BoxWalls") = {1,3,4,2};
Physical Curve("Ellipse") = {6, 5, 8, 7};
Physical Surface("Fluid") = {2};
Physical Surface("Ell") = {1};

Rotate {{0, 0, 1}, {0.5, 0.08, 0}, Pi/4} {
  Surface{1}; 
}