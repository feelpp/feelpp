h=0.1;
a=-0.5;
b=1;
c=0;
d=2;

Point(1) = {a, c, 0., h};
Point(2) = {b, c, 0., h};
Point(3) = {b, d, 0., h};
Point(4) = {a, d, 0., h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};

Physical Line("top") = {3};
Physical Line("bottom") = {1};
Physical Line("left") = {4};
Physical Line("right") = {2};
Physical Surface("Omega") = {1};
