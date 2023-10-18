// Gmsh project created on Fri Dec 16 10:05:59 2022

h=0.1;

SetFactory("OpenCASCADE");

d=0.25; 
L=2;
D=0.5;

L_ext=2.5;
D_ext=1;

//+
Cylinder(1) = {0, 0, 0, 0, 0, 1, D, 2*Pi};
//+
Cylinder(2) = {0, 0, 0, 0, 0, 1, D_ext, 2*Pi};
//+
Dilate {{0, 0, 0}, {1, 1, L}} {
  Volume{1}; 
}
//+
Dilate {{0, 0, 0}, {1, 1, L_ext}} {
  Volume{2}; 
}
//+
Translate {0, 0, 0.4} {
  Volume{1}; 
}
//+
Surface Loop(6) = {10, 11, 12};
//+
Surface Loop(7) = {13, 14, 15};
//+
Circle(16) = {0, 0, 2.5, d, 0, 2*Pi};

//+
Extrude {0, 0, -0.1} {
  Curve{16}; 
}
//+
Curve Loop(17) = {13};
//+
Curve Loop(18) = {18};
//+
Plane Surface(17) = {17, 18};
//+
Surface Loop(8) = {11, 10, 12};
//+
Curve Loop(20) = {10};
//+
Curve Loop(21) = {16};
//+
Plane Surface(18) = {20, 21};
//+
Surface Loop(9) = {18, 16, 17, 10, 12, 13, 15};
//+
Volume(3) = {9};

Delete{Volume{1:2};}

Characteristic Length {PointsOf{Surface{13};Surface{17};Surface{15};}} = 0.5*h;

Characteristic Length {PointsOf{Surface{10};Surface{12};Surface{18};}} = 2*h;
//+
Physical Volume("Insulation") = {3};
//+
Physical Surface("CavitySide") = {13};
//+
Physical Surface("CavityBottom") = {15};
//+
Physical Surface("CavityTop") = {17};
