Mesh.RecombineAll = 1;

Point(1) = {1.6666666666666667, 0, -0.1, 1};
Point(2) = {0, 1.6666666666666667, -0.1, 1};
Point(3) = {-1.6666666666666667, 0, -0.1, 1};
Point(4) = {0, -1.6666666666666667, -0.1, 1};
Point(5) = {5.0, 0, -0.1, 1};
Point(6) = {0, 5.0, -0.1, 1};
Point(7) = {-5.0, 0, -0.1, 1};
Point(8) = {0, -5.0, -0.1, 1};
Point(9) = {0, 0, -0.1, 1};

Circle(1) = {5, 9, 6};
Circle(2) = {6, 9, 7};
Circle(3) = {7, 9, 8};
Circle(4) = {8, 9, 5};
Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9) = {1, 2};
Line(10) = {2, 3};
Line(11) = {3, 4};
Line(12) = {4, 1};

Line Loop(1) = {9, 10, 11, 12};
Line Loop(2) = {1, -6, -9, 5};
Line Loop(3) = {2, -7, -10, 6};
Line Loop(4) = {3, -8, -11, 7};
Line Loop(5) = {4, -5, -12, 8};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

Transfinite Line {1, 2, 3, 4, 9, 10, 11, 12} = 7;
Transfinite Line {5, 6, 7, 8} = 5;
Transfinite Surface {1} = {1, 2, 3, 4};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};
Recombine Surface {1, 2, 3, 4, 5};

center[] = Extrude {0, 0, 0.2} {
  Surface{1};
  Layers{1};
  Recombine;
};
s1[] = Extrude {0, 0, 0.2} {
  Surface{2};
  Layers{1};
  Recombine;
};
s2[] = Extrude {0, 0, 0.2} {
  Surface{3};
  Layers{1};
  Recombine;
};
s3[] = Extrude {0, 0, 0.2} {
  Surface{4};
  Layers{1};
  Recombine;
};
s4[] = Extrude {0, 0, 0.2} {
  Surface{5};
  Layers{1};
  Recombine;
};

Physical Surface("ZMoins") = {1, 2, 3, 4, 5};
Physical Surface("ZPlus") = {center[0], s1[0], s2[0], s3[0], s4[0]};
Physical Surface("Outer") = {s1[2], s2[2], s3[2], s4[2]};
Physical Volume("Shell") = {center[1], s1[1], s2[1], s3[1], s4[1]};
