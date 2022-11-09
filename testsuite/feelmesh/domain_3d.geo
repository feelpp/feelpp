// Gmsh project created on Wed May 11 18:17:21 2022
SetFactory("OpenCASCADE");
//+
h = 0.05;
Mesh.CharacteristicLengthMax = h;
//+
Box(1) = {0, 0, 0, 1, 1, 1};
Box(2) = {0,0,0,0.5,0.5,0.5};
b() = BooleanFragments{ Volume{1}; Delete; }{Volume{2}; Delete;};
Physical Volume("MatOne") = {3};
Physical Surface("RequiredBoundaryOfRequiredElements") = {6,7,9};
Physical Surface("OtherRequiredBoundary") = {8};