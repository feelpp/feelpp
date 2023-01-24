L = 0.4;
A = L/4; // large axis
B = L/8; // small axis
lc=0.01;
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {L, 0, 0, lc};
//+
Point(3) = {L, 7*L, 0, lc};
//+
Point(4) = {0, 7*L, 0, lc};
//+
Point(5) = {0.5*L, 6*L, 0, lc};
//+
Point(6) = {0.5*L-A*0.5, 6*L, 0, lc/2};
//+
Point(7) = {0.5*L+A*0.5, 6*L, 0, lc/2};
//+
Point(8) = {0.5*L, 6*L+0.5*B, 0, lc/2};
//+
Point(9) = {0.5*L, 6*L-0.5*B, 0, lc/2};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Ellipse(5) = {7, 5, 8};
//+
Ellipse(6) = {8, 5, 6};
//+
Ellipse(7) = {6, 5, 9};
//+
Ellipse(8) = {9, 5, 7};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(1) = {2};
//+
Plane Surface(2) = {1, 2};
//+
Physical Curve("Top") = {3};
//+
Physical Curve("Walls") = {4, 2};
//+
Physical Curve("EllipseSurface") = {6, 5, 8, 7};
//+
Physical Curve("Bottom") = {1};
//+
Physical Surface("Fluid") = {2};
//+
Physical Surface("EllipseVolume") = {1};

//+
Rotate {{0, 0, 1}, {0.2, 2.4, 0}, Pi/4} {
  Curve{6}; Point{8}; Point{6}; Curve{7}; Point{9}; Curve{8}; Point{7}; Curve{5}; 
}
