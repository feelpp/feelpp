// Gmsh project created on Wed May 11 18:17:21 2022
SetFactory("OpenCASCADE");
//+
h = 0.05;
Mesh.CharacteristicLengthMax = h;
//+
Sphere(1) = {0, 0, 0, 1,-Pi/2,Pi/2,2*Pi};
Sphere(2) = {0,0,0,0.25,-Pi/2,Pi/2,2*Pi};
b() = BooleanFragments{ Volume{1}; Delete; }{Volume{2}; Delete;};
Physical Volume("MatTwo") = {2};
Physical Volume("MatOne") = {3};
Physical Surface("RequiredBoundaryOfRequiredElements") = {2};
Physical Surface("OtherRequiredBoundary") = {3};