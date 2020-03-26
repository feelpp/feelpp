h = 0.1;
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
Point(1) = {xmin,ymin,0.0,h};
Point(2) = {xmax,ymin+0,0.0,h};
Point(3) = {xmax+0,ymax,0.0,h};
Point(4) = {xmin+0,ymax+0,0.0,h};
Point(5) = {xmin, ymin +0.2,0.0,h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Line Loop(6) = {4, 5, 1, 2, 3};
Plane Surface(7) = {6};

Physical Line("Gamma0") = {2,3,4};
Physical Line("Gamma1") = {1,5};
Physical Surface("Omega") = {7};
