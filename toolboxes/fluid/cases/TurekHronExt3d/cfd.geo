h=0.03; // M1:0.03
myhA=h;
myhB=h/2.;
myhC=h/4.;
Long=2.5;
Haut=0.41;
Circle_xc=0.2;
Circle_yc=0.2;
lgstruct=0.35101;
hstruct=0.02;
Circle_radius=0.05;
//Circle_r=;
Point(1) = {0., 0., 0., myhA};
Point(2) = {Long, 0., 0., myhA};
Point(3) = {Long, Haut, 0., myhA};
Point(4) = {0., Haut, 0., myhA};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Point(5) = {Circle_xc, Circle_yc, 0., myhC};
Point(6) = {0.6-lgstruct,0.19,0.,myhC};
Point(7) = {0.6-lgstruct,0.19+hstruct,0.,myhC};
Point(8) = {Circle_xc-Circle_radius,Circle_yc,0.,myhC};
Circle(5) = {7, 5, 8};
Circle(6) = {8, 5, 6};
Point(9) = {0.6,0.19,0.,myhC};
Point(10) = {0.6,0.19+hstruct,0.,myhC};
Point(11) = {0.6,0.19+hstruct/2.,0,myhC};
Point(12) = {1.5,0.19+hstruct/2.,0,myhC};
Point(13) = {1.5,0.19+ Haut/4.,0,myhB};
Point(14) = {1.5,0.19- Haut/4.,0,myhB};
Point(15) = {0.075,0.19+ Haut/4.,0,myhB};
Point(16) = {0.075,0.19- Haut/4.,0,myhB};
Point(17) = {0.075,0.19+hstruct/2.,0,myhC};

Line(7) = {6, 9};
Line(8) = {9, 11};
Line(9) = {11, 10};
Line(10) = {10, 7};
Line(11) = {15, 13};
Line(12) = {13, 12};
Line(13) = {12, 14};
Line(14) = {14, 16};
Line(15) = {16, 17};
Line(16) = {17, 15};
Line(17) = {17, 8};
Line(18) = {11, 12};


Line Loop(11) = {3, 4, 1, 2};
Line Loop(13) = {11, 12, 13, 14, 15, 16};
Plane Surface(1) = {11, 13};
Line Loop(14) = {11, 12, -18, 9, 10, 5, -17, 16};
Plane Surface(2) = {14};
Line Loop(15) = {15, 17, 6, 7, 8, 18, 13, 14};
Plane Surface(3) = {15};

Extrude {0, 0, Haut} {
  Surface{1,2,3};
}

Physical Surface("inlet") = {37};
Physical Surface("outlet") = {45};
Physical Surface("wall2") = {103, 133, 99, 137, 95, 141};
Physical Surface("wall1") = {33, 70, 112, 154, 41, 1, 3, 2};
Physical Volume("Fluid") = {3, 2, 1};
