h = 1;

xmin = 0;
ymin = 0;
zmin = 0;

xmax = 1;
ymax = 1;
zmax = 1;
Point(1) = {xmin,ymin,zmin,h};
Point(2) = {xmax,ymin,zmin,h};
Point(3) = {xmax,ymax,zmin,h};
Point(4) = {xmin,ymax,zmin,h};
Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,1};
Line(4) = {1,2};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
Transfinite Surface{6} = {1,2,3,4};

out[] = Extrude {0,0,zmax-zmin} {
  Surface{6};
  Layers { {(zmax-zmin)/h}, {1.0} };
};

Physical Surface("Border") = {6, 15, 19, 23, 27, 28};
Physical Volume("Cube") = {1};
