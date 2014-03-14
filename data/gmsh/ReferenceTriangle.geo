// -*- Mode: c++ -*-
h=4;
Point(1) = {-1,-1,0,h};
Point(2) = {1,-1,0,h};
Point(3) = {-1,1,0,h};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};
Line Loop(4) = {3,1,2};
Plane Surface(5) = {4};
Physical Line("hor") = {1};
Physical Line("hypo") = {2};
Physical Line("vert") = {3};
Physical Surface(9) = {5};
