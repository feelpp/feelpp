h = 0.019;
xmin = -1;
xmax =  0.5;
ymin = -1;
ymax =  0.5;
zmin = -1;
zmax =  1;
Point(1) = {xmin,ymin,zmin,h};
Point(2) = {xmax,ymin,zmin,h};
Line(1) = {1,2};

Extrude Line {1, {0,ymax-ymin,0} } {
  Layers { (ymax-ymin)/h};
};

Extrude Surface {5, {0,0,zmax-zmin} } {
  Layers { (zmax-zmin)/h};
};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

Physical Line( 7) = { 7};
Physical Line( 8) = { 8};
Physical Line( 9) = { 9};
Physical Line(10) = {10};

Physical Line(12) = {12};
Physical Line(13) = {13};
Physical Line(17) = {17};
Physical Line(21) = {21};
Physical Surface("Border") = { 5,22,27,14,26,18};
Physical Volume("firstMat") = {1};
