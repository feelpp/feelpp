h=0.01; //change by user

Long=2.5;
Haut=0.41;
Circle_xc=0.2;
Circle_yc=0.2;
lgstruct=0.35101;
//lgstruct=0.35;
hstruct=0.02;
Circle_r=0.05;

Point(1) = {Circle_xc, Circle_yc, 0.,h};
Point(2) = {0.6-lgstruct,0.19,0.,h};
Point(3) = {0.6-lgstruct,0.19+hstruct,0.,h};

Point(4) = {0.6,0.19,0.,h};
Point(5) = {0.6,0.19+hstruct,0.,h};

Line(1) = {2, 4};
Line(2) = {4, 5};
Line(3) = {5, 3};
Circle(4) = {2, 1, 3};
Line Loop(5) = {3, -4, 1, 2};
Plane Surface(1) = {5};

Physical Line("fixed-wall") = {4};
Physical Line("free-wall") = {1,2,3};
Physical Surface("beam") = {1};
