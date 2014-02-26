h = 0.005;
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
Point(1) = {xmin,ymin,0.0,h};
Point(2) = {xmax,ymin+0,0.0,h};
Point(3) = {xmax+0,ymax,0.0,h};
Point(4) = {xmin+0,ymax+0,0.0,h};
Line(1) = {4,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
tmp[] = Extrude {0,0.0,1} {
  Surface{6};
};
Physical Volume(1) = tmp[1];

Physical Surface(6) = {6};
Physical Surface("wall2") = {27};
Physical Surface("wall") = {15, 28, 6, 23, 19};

