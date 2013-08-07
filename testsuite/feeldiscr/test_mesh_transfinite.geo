h = 10;

// Mesh.CharacteristicLengthExtendFromBoundary=1;
// Mesh.CharacteristicLengthFromPoints=1;

Point(1) = {0, 0, 0, h};
Point(2) = {10, 0, 0, h};
Point(3) = {10, 10, 0, h};
Point(4) = {0, 10, 0, h};
Point(5) = {0, 0, 10, h};
Point(6) = {10, 0, 10, h};
Point(7) = {10, 10, 10, h};
Point(8) = {0, 10, 10, h};
Line(1) = {4, 8};
Line(2) = {8, 7};
Line(3) = {7, 3};
Line(4) = {3, 4};
Line(5) = {8, 5};
Line(6) = {4, 1};
Line(7) = {3, 2};
Line(8) = {7, 6};
Line(9) = {6, 2};
Line(10) = {2, 1};
Line(11) = {1, 5};
Line(12) = {5, 6};
Line Loop(14) = {1, 2, 3, 4};
Plane Surface(14) = {14};
Line Loop(16) = {6, 11, -5, -1};
Plane Surface(16) = {16};
Line Loop(18) = {10, -6, -4, 7};
Plane Surface(18) = {18};
Line Loop(20) = {3, 7, -9, -8};
Plane Surface(20) = {20};
Line Loop(22) = {2, 8, -12, -5};
Plane Surface(22) = {22};
Line Loop(24) = {10, 11, 12, 9};
Plane Surface(24) = {24};
Surface Loop(26) = {16, 18, 24, 22, 14, 20};
Volume(26) = {26};


Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12} = 2 Using Progression 1;
Transfinite Surface {14,16,18,20,22,24};
Transfinite Volume {26};

Physical Surface("Bottom") = {24};
Physical Surface("Right") = {20};
Physical Surface("Top") = {14};
Physical Surface("Left") = {16};
Physical Surface("Front") = {18};
Physical Surface("Back") = {22};

Physical Volume("OmegaFluide") = {26};
