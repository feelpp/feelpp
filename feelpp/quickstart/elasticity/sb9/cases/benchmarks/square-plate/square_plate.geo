h = 1.25;
h = DefineNumber[ h, Name "Parameters/h" ];
L = 10.0;
W = 10.0;
T = 0.2;
nx = Ceil( L/h );
ny = Ceil( W/h );
nz = 1;
If ( nx < 8 )
  nx = 8;
EndIf
If ( ny < 8 )
  ny = 8;
EndIf

Mesh.RecombineAll = 1;

Point(1) = {0, 0, -T/2, h};
Point(2) = {L, 0, -T/2, h};
Point(3) = {0, W, -T/2, h};
Point(4) = {L, W, -T/2, h};
Point(5) = {0, 0, T/2, h};
Point(6) = {L, 0, T/2, h};
Point(7) = {0, W, T/2, h};
Point(8) = {L, W, T/2, h};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line(5) = {5, 6};
Line(6) = {6, 8};
Line(7) = {8, 7};
Line(8) = {7, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {4, 8};
Line(12) = {3, 7};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Line Loop(3) = {1, 10, -5, -9};
Plane Surface(3) = {3};
Line Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};
Line Loop(5) = {-3, 11, 7, -12};
Plane Surface(5) = {5};
Line Loop(6) = {-4, 12, 8, -9};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Line {1, 3, 5, 7} = nx + 1;
Transfinite Line {2, 4, 6, 8} = ny + 1;
Transfinite Line {9, 10, 11, 12} = nz + 1;
Transfinite Surface {1} = {1, 2, 4, 3};
Transfinite Surface {2} = {5, 6, 8, 7};
Transfinite Surface {3} = {1, 2, 6, 5};
Transfinite Surface {4} = {2, 4, 8, 6};
Transfinite Surface {5} = {3, 4, 8, 7};
Transfinite Surface {6} = {1, 3, 7, 5};
Transfinite Volume {1} = {1, 2, 4, 3, 5, 6, 8, 7};

Recombine Surface {1, 2, 3, 4, 5, 6};
Recombine Volume {1};

Physical Surface("ZMoins") = {1};
Physical Surface("ZPlus") = {2};
Physical Surface("YMoins") = {3};
Physical Surface("XPlus") = {4};
Physical Surface("YPlus") = {5};
Physical Surface("XMoins") = {6};
Physical Volume("Shell") = {1};

// P1..P8 follow the MATLAB/PDF numbering used in validation.pdf.
Physical Point("P1") = {1};
Physical Point("P2") = {2};
Physical Point("P3") = {3};
Physical Point("P4") = {4};
Physical Point("P5") = {5};
Physical Point("P6") = {6};
Physical Point("P7") = {7};
Physical Point("P8") = {8};
