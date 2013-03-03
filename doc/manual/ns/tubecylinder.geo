h = 0.05;
Mesh.RecombinationAlgorithm=0;
xi = 0;
xf = 2.2;
yi =0;
yf = 0.41;

xc = 0.2;
yc = 0.2;
R = 0.05;

Point(1) = {xi,yi,0.0,h};
Point(2) = {xf,yi,0.0,h};
Point(3) = {xf,yf,0.0,h};
Point(4) = {xi,yf,0.0,h};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(17) = {1,2,3,4};

Point(5) = {xc,yc,0.0,h/6};
Point(6) = {xc-R,yc,0.0,h/6};
Point(7) = {xc,yc+R,0.0,h/6};
Point(8) = {xc+R,yc,0.0,h/6};
Point(9) = {xc,yc-R,0.0,h/6};
Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};
Line Loop(18) = {5,6,7,8};


Plane Surface(18) = {17,18};
Physical Line("inlet") = {4};
Physical Line("outlet") = {2};
Physical Line("wall") = {1,3};
Physical Line("cylinder") = {5,6,7,8};
Physical Surface(19) = {18};
