h=1;

x0 = -30;
y0 =   0;
z0 = -10;

x1 = 10;
y1 = 10;
z1 =  0;

p1 = newp; Point(p1) = {x0, y0, z1, h};
p2 = newp; Point(p2) = {x1, y0, z1, h};
p3 = newp; Point(p3) = {x1, y1, z1, h};
p4 = newp; Point(p4) = {x0, y1, z1, h};

l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};
ll1 = newll; Line Loop(ll1) = {l1:l4};

p5 = newp; Point(p5) = {-0.25, 1.5,  z1, h};
p6 = newp; Point(p6) = { 0.25, 1.5,  z1, h};
p7 = newp; Point(p7) = { 0.25, 2.0,  z1, h};
p8 = newp; Point(p8) = {-0.25, 2.0,  z1, h};
l5 = newl; Line(l5) = {p5, p6};
l6 = newl; Line(l6) = {p6, p7};
l7 = newl; Line(l7) = {p7, p8};
l8 = newl; Line(l8) = {p8, p5};
ll2 = newll; Line Loop(ll2) = {l5:l8};

Plane Surface(1) = {ll1, -ll2};
Plane Surface(2) = {ll2};

Translate {0, 0, -0.1} { Duplicata{ Surface{2}; } }
Translate {0, 0, -(z1-z0)} { Duplicata{ Surface{1}; } }

l9 = newl; Line(l9) = {19, 1};
l10 = newl; Line(l10) = {28, 4};
l11 = newl; Line(l11) = {20, 2};
l12 = newl; Line(l12) = {24, 3};

ll3 = newll; Line Loop(ll3) = {l4, -l9, -20, l10};
Plane Surface(3) = {ll3};
ll4 = newll; Line Loop(ll4) = {-l2, -l11, 18, l12};
Plane Surface(4) = {ll4};
ll5 = newll; Line Loop(ll5) = {l1, -l11, -17, l9};
Plane Surface(5) = {ll5};
ll6 = newll; Line Loop(ll6) = {-l3, -l12, 19, l10};
Plane Surface(6) = {ll6};


l13 = newl; Line(l13) = {5, 9};
l14 = newl; Line(l14) = {8, 18};
l15 = newl; Line(l15) = {6, 10};
l16 = newl; Line(l16) = {7, 14};

ll7 = newll; Line Loop(ll7) = {12, -35, -6, 33};
Plane Surface(7) = {ll7};
ll8 = newll; Line Loop(ll8) = {-35,7, 36, -13};
Plane Surface(8) = {ll8};
ll9 = newll; Line Loop(ll9) = {-14,-36,8,34};
Plane Surface(9) = {ll9};
ll10 = newll; Line Loop(ll10) = {15,-33,-9,34};
Plane Surface(10) = {ll10};

ll11 = newll; Line Loop(ll11) = {21,22,23,24};
Plane Surface(12) = {ll11};

Surface Loop(41) = {2, 8, 7, 11, 9, 10};
Volume(42) = {41};

Surface Loop(42) = {1,3,5,4,6,16,7,8,9,10,11,12};
Volume(43) = {42};

Physical Surface("Dirichlet") = {4,3};
Physical Volume("v1") = {42};
Physical Volume("v2") = {43};
