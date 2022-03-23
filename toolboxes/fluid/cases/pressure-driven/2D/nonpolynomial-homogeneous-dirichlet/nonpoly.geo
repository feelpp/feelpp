h=0.1;
xstart=0;
xlength=1;
ystart=0;
yheight=1;
Point(1) = {xstart, ystart, 0., h};
Point(2) = {xstart+xlength, ystart, 0., h};
Point(3) = {xstart+xlength, ystart+yheight, 0., h};
Point(4) = {xstart, ystart+yheight, 0., h};
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
