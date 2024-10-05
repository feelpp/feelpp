// -*- mode: c¯ -*-
h=0.1;

L=1;
H=1;
Point(1)={0, 0, 0, h};
Point(2)={L, 0, 0, h};
Point(3)={L, H, 0, h};
Point(4)={0, H, 0, h};

Point(5)={0.5,1,0,h};
Point(6)={0,0.5,0,h};
Point(7)={1,0.5,0,h};
Point(8)={0.5,0,0,h};
Point(9)={0.5,0.5,0,h};
Line(1)={1,8};
Line(2)={8,2};
Line(3)={2,7};
Line(4)={7,3};
Line(5)={3,5};
Line(6)={5,4};
Line(7)={4,6};
Line(8)={6,1};
Line(9)={8,9};
Line(10)={9,5};
Line(11)={6,9};
Line(12)={9,7};
Line Loop(14)={7,11,10,6};
Plane Surface(14)={14};
Line Loop(16)={8,1,9,-11};
Plane Surface(16)={16};
Line Loop(18)={2,3,-12,-9};
Plane Surface(18)={18};
Line Loop(20)={12,4,5,-10};
Plane Surface(20)={20};

//physicalentities
Physical Line("Tflux")={1, 2};
Physical Line("Tfourier")={3,4,5,6,7,8};
Physical Surface("Fin0")={16};
Physical Surface("Fin1")={18};
Physical Surface("Fin2")={14};
Physical Surface("Fin3")={20};
