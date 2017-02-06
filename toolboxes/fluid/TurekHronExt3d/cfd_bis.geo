h=0.1;//0.03;
myhext=h;
//myhprecis=0.003;
myhprecis=myhext/5.;
myhprecis=myhext/3.;
Long=2.5;
Haut=0.41;
Circle_xc=0.2;
Circle_yc=0.2;
lgstruct=0.35101;
hstruct=0.02;
Circle_radius=0.05;
//Circle_r=;
Point(1) = {0., 0., 0., myhext};
Point(2) = {Long, 0., 0., myhext};
Point(3) = {Long, Haut, 0., myhext};
Point(4) = {0., Haut, 0., myhext};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Point(5) = {Circle_xc, Circle_yc, 0., myhprecis};

Point(6) = {Circle_xc+Circle_radius,Circle_yc,0.,myhprecis};
Point(8) = {Circle_xc-Circle_radius,Circle_yc,0.,myhprecis};
Point(9) = {Circle_xc,Circle_yc-Circle_radius,0.,myhprecis};
Point(10) = {Circle_xc,Circle_yc+Circle_radius,0.,myhprecis};
Circle(5) = {6, 5, 10};
Circle(6) = {10, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};


Line Loop(9) = {3, 4, 1, 2};
Line Loop(10) = {5, 6, 7, 8};
Plane Surface(11) = {9, 10};
Extrude {0, 0, Haut} {
  Surface{11};
}

Physical Surface("inlet") = {28};
Physical Surface("outlet") = {36};
Physical Surface("wall1") = {24, 11, 53, 32};
Physical Surface("wall2") = {52, 40, 44, 48};
Physical Volume("Fluid") = {1};

