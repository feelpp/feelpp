// Gmsh project created on Sat Dec 14 17:02:12 2013

h=0.01;

Point(1) = {0,0,0,h};
Point(2) = {0,-0.05,0,h};
Point(3) = {0.05,0,0,h};
Point(4) = {0,0.05,0,h};
Point(5) = {-0.05,0,0,h};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Extrude {0,0,0.5} {Surface{6}; Layers{40};}

Physical Surface("inlet") = {6};
Physical Surface("outlet") = {28};
Physical Surface("wall") = {15,19,23,27};
Physical Volume(29) = {1};
