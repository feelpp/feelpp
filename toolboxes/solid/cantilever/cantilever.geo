h = 0.01;
a = 0.1;
h = DefineNumber[ h, Name "Parameters/h" ];
a = DefineNumber[ a, Name "Parameters/a" ];
Point(1) = {0, -a/2, 0, h};
Point(2) = {0, a/2, 0, h};
Line(1) = {2, 1};
Extrude {0, 0, 2*a} {
  Line{1};
}
Extrude {10*a, 0, 0} {
  Surface{5};
}
Physical Surface("fixed") = {5};
Physical Surface("load") = {27};
Physical Volume("cantilever") = {1};

