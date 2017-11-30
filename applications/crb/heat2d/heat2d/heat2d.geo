// -*- mode: c++ -*-
/// [geo]
h=0.1;
a=0;
b=1;
c=2;
d=1;
Point(1) = {a,0,0,h}; 
Point(2) = {b,0,0,h}; 
Point(3) = {c,0,0,h}; 
Point(4) = {c,d,0,h}; 
Point(5) = {b,d,0,h}; 
Point(6) = {a,d,0,h}; 


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(8) = {2, 5};
Line Loop(5) = {1,8,5,6};
Line Loop(6) = {2,3,4,-8};

Plane Surface(5)={5};
Plane Surface(6)={6};

Physical Line("BR") = {3};
Physical Line("BL") = {6};

Physical Surface("SL") = {5};
Physical Surface("SR") = {6};
/// [geo]

//Physical Surface("Omega") = {5,6};
