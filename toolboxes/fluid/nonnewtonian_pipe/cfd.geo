h=1;//0.05*1e-3;
Long=1e-3;
Haut=4e-3;
Point(1) = {0., 0., 0., h};
Point(2) = {Long, 0., 0., h};
Point(3) = {Long, Haut, 0., h};
Point(4) = {0., Haut, 0., h};
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
