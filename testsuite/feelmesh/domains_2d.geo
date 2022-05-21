// Gmsh project created on Wed May 11 18:17:21 2022
SetFactory("OpenCASCADE");
//+
h = 0.05;
Mesh.CharacteristicLengthMax = h;
//+
Rectangle(1) = {0, 0, 0, 0.75, 0.75, 0.};
Rectangle(2) = {0,0,0,0.5,0.5,0.};
b() = BooleanFragments{ Surface{1}; Delete; }{Surface{2}; Delete;};
Physical Surface("MatTwo") = {2};
Physical Surface("MatOne") = {3};
Physical Curve("RequiredBoundaryOfRequiredElements") = {1,2};
Physical Curve("OtherRequiredBoundary") = {4};