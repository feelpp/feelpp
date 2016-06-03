h=0.25;

x0 = 0;
y0 = 0;
z0 = 0;

x1 = 1;
y1 = 1;
z1 = 1;


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

out[] = Extrude{0,0,z1-z0}{ Surface{1}; };
//Physical Surface("OxOy") = {1};
//Physical Surface("OxOy+") = {26};
//Physical Surface("OyOz") = {25};
//Physical Surface("OyOz+") = {17};
//Physical Surface("OxOz") = {13};
//Physical Surface("OxOz+") = {21};
Physical Surface("Border") = {1,26,25,17,13,21};
Physical Volume("COIL") = {1};
