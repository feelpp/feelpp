h = 0.5;
h = DefineNumber[ h, Name "Parameters/h" ];

L = 2;
W = 1;
T = 0.1;

Mesh.RecombineAll = 1;

Point(1) = {0, 0, -T/2, h};
Point(2) = {L, 0, -T/2, h};
Point(3) = {L, W, -T/2, h};
Point(4) = {0, W, -T/2, h};
Point(5) = {0, 0,  T/2, h};
Point(6) = {L, 0,  T/2, h};
Point(7) = {L, W,  T/2, h};
Point(8) = {0, W,  T/2, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Curve Loop(3) = {1, 10, -5, -9};
Plane Surface(3) = {3};
Curve Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};
Curve Loop(5) = {3, 12, -7, -11};
Plane Surface(5) = {5};
Curve Loop(6) = {4, 9, -8, -12};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Line {1, 3, 5, 7} = 5;
Transfinite Line {2, 4, 6, 8} = 4;
Transfinite Line {9, 10, 11, 12} = 2;
Transfinite Surface {1} = {1, 2, 3, 4};
Transfinite Surface {2} = {5, 6, 7, 8};
Transfinite Surface {3} = {1, 2, 6, 5};
Transfinite Surface {4} = {2, 3, 7, 6};
Transfinite Surface {5} = {3, 4, 8, 7};
Transfinite Surface {6} = {4, 1, 5, 8};
Transfinite Volume {1} = {1, 2, 3, 4, 5, 6, 7, 8};

Recombine Surface {1, 2, 3, 4, 5, 6};
Recombine Volume {1};

Physical Surface("Dirichlet") = {1, 2, 3, 4, 5, 6};
Physical Volume("Shell") = {1};
