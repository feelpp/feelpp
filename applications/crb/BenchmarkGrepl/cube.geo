h=0.1;

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

out[]  = Extrude{0,0,-(z1-z0)}{ Surface{1}; };
out[] += Extrude{0,0,-0.1}{ Surface{2}; };
out[] += Extrude{0,0,-(z1-z0)+0.1}{ Surface{74}; };

sl1 = newsl; Surface Loop(sl1) = {1,23,52,31,35,27,96,61,65,69,73,74};
v = newv; Volume(v) = {sl1};

Physical Volume("v1") = {v};
Physical Volume("v2") = {out[11]};


