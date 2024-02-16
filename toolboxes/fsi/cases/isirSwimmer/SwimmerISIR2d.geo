h=0.02;
myhext=h;
//myhprecis=0.001;
myhprecis=myhext/200;
Long=0.2;
Haut=0.2;
Circle_xc=0.1;
Circle_yc=0.1;
lgstruct=0.008; // swimmer length
hstruct=0.0005; // magnet height
Circle_radius=0.00075; //magnet radius
Tail_radius=0.00008;
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
Point(8) = {Circle_xc-hstruct,Circle_yc+Circle_radius,0.,myhprecis};
Point(9) = {Circle_xc-hstruct,Circle_yc-Circle_radius,0.,myhprecis};
Point(10) = {Circle_xc-hstruct/2,Circle_yc,0.,myhprecis};
Point(11) = {Circle_xc+lgstruct,Circle_yc,0,myhprecis};


// surface fluid
Line Loop(11) = {3, 4, 1, 2};
Line Loop(12) = {5, 6, 7, 8, 9, 10,11,-11};

// surface structure

Line Loop(15) = {10, 5, 6, 7, 8, 9};

//Physical Line("FEELPP_GMSH_PHYSICALNAME_IGNORED") = {11};

Physical Line("fluid-inlet") = {4};
Physical Line("fluid-wall") = {1,3};
Physical Line("fluid-outlet") = {2};

Point(16) = {Circle_xc, Circle_yc+Circle_radius, 0, myhprecis};
Point(17) = {Circle_xc, Circle_yc-Circle_radius, 0, myhprecis};

Line(29) = {16, 8};
Line(30) = {8, 9};
Line(31) = {9, 17};
Line(32) = {17, 16};
Line Loop(33) = {30, 31, 32, 29};
Plane Surface(34) = {33};
Point(22) = {Circle_xc+lgstruct-Tail_radius, 0.1, 0, myhprecis};
Point(23) = {Circle_xc+lgstruct-Tail_radius, 0.1+Tail_radius, 0, myhprecis};
Point(24) = {Circle_xc+lgstruct-Tail_radius, 0.1-Tail_radius, 0, myhprecis};
Circle(35) = {23, 22, 11};
Circle(36) = {11, 22, 24};

Line(37) = {16, 23};

Line(38) = {24, 17};
Line Loop(39) = {37, 35, 36, 38, 32};
Plane Surface(40) = {39};
Line Loop(41) = {38, -31, -30, -29, 37, 35, 36};
Plane Surface(42) = {11, 41};
Physical Line("Tail-boundary") = {38, 37, 35, 36};
Physical Line("Magnet-boundary") = {31, 30, 29};
Physical Line("Tail-Magnet-interface") = {32};
Physical Surface("Fluid") = {42};
Physical Surface("Magnet") = {34};
Physical Surface("Elastic") = {40};
