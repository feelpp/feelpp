// Gmsh project created on Fri Jun 16 11:48:55 2017
SetFactory("OpenCASCADE");

//+
Cylinder(1) = {0, 0, 0, 0, 0, 4.61, 53.2, 1.99*Pi};
//+
Cylinder(2) = {0, 0, 0, 0, 0, 4.61, 30.6, 1.99*Pi};

BooleanDifference(3) = {Volume{1}; Delete;} {Volume{2}; Delete;};
//+
//+
Sphere(4) = {0, 0, 0, 20, -Pi/2, Pi/2, 2*Pi};
//+
// Sphere(5) = {0, 0, 0, 200, -Pi/2, Pi/2, 2*Pi};
//+
// Sphere(6) = {0, 0, 0, 300, -Pi/2, Pi/2, 2*Pi};

// BooleanDifference(7) = {Volume{6}; Delete;} {Volume{5};};
// BooleanDifference(8) = {Volume{5}; Delete;} {Volume{3,4};};

// Physical Volume("inf") = {7};
// Physical Volume("ext") = {8};
Physical Volume("cond") = {3};
Physical Volume("mgn") = {4};
Physical Surface("in") = {3};
Physical Surface("out") = {5};
Physical Surface("Rint") = {6};
Physical Surface("Rext") = {1};
Physical Surface("top") = {2};
Physical Surface("bottom") = {4};
// Physical Surface("Sext") = {8};
// Physical Surface("Sinf") = {9};

h=0.6;
Characteristic Length {1,2,3,4,5,6,7,8} = 5*h;
Characteristic Length {9,10} = 3*h;
// Characteristic Length {11,12} = 10*h;
// Characteristic Length {13,14} = 20*h;

