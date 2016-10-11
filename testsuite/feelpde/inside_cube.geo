h=0.25;
eps=0.001;

x0 = 0;
y0 = 0;
z0 = 0;

x1 = 1;
y1 = 1;
z1 = 0.5-eps;
z2 = 0.5+eps;
z3 = 1;

xx0 = 0.2;
yy0 = 0.2;
xx1 = 0.8;
yy1 = 0.8;

Point(1) = {x0, y0, z0, h};
Point(2) = {x1, y0, z0, h};
Point(3) = {x1, y1, z0, h};
Point(4) = {x0, y1, z0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1:4};
Plane Surface(1) = {1};

Point(10) = {xx0, yy0, z1, h};
Point(20) = {xx1, yy0, z1, h};
Point(30) = {xx1, yy1, z1, h};
Point(40) = {xx0, yy1, z1, h};

Line(10) = {10, 20};
Line(20) = {20, 30};
Line(30) = {30, 40};
Line(40) = {40, 10};
Line Loop(10) = {10,20,30,40};
Plane Surface(2) = {10};

out1[] = Extrude{0,0,z3-z0}{ Surface{1}; };
out2[] = Extrude{0,0,z2-z1}{ Surface{2}; };

Physical Surface("Border") = {1,out1[0],out1[2],out1[3],out1[4],out1[5]};
Physical Volume("COIL") = {out1[1],out2[1]};
