h=0.03;
myhext=h;
//myhprecis=0.003;
myhprecis=myhext/5.;
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
Point(6) = {0.6-lgstruct,0.19,0.,myhprecis};
Point(7) = {0.6-lgstruct,0.19+hstruct,0.,myhprecis};
Point(8) = {Circle_xc-Circle_radius,Circle_yc,0.,myhprecis};
Circle(5) = {7, 5, 8};
Circle(6) = {8, 5, 6};
Point(9) = {0.6,0.19,0.,myhprecis};
Point(10) = {0.6,0.19+hstruct,0.,myhprecis};

Point(11) = {0.6,0.19+hstruct/2.,0,myhprecis};
Point(12) = {1.5,0.19+hstruct/2.,0,myhprecis*3.};

Line(7) = {6, 9};
Line(8) = {9, 11};
Line(9) = {11, 10};
Line(10) = {10, 7};
Line(11) = {11, 12};


// surface fluid
Line Loop(11) = {3, 4, 1, 2};
Line Loop(12) = {5, 6, 7, 8, 9, 10,11,-11};
Plane Surface(1) = {11, 12};

// surface structure
Circle(14) = {7, 5, 6};
Line Loop(15) = {10, 14, 7, 8, 9};
Plane Surface(2) = {15};

Line Loop(16) = {5, 6, -14};
Plane Surface(17) = {16};

Physical Line("line-downstream") = {11};

Physical Line("inlet") = {4};
Physical Line("wall") = {1,3};
Physical Line("fluid-cylinder") = {5,6};
Physical Line("outlet") = {2};
Physical Line("fluid-beam") = {7,8,9,10};
Physical Line("cylinder-beam") = {14};

Physical Surface("Fluid") = {1};
Physical Surface("beam") = {2};
Physical Surface("cylinder") = {17};

