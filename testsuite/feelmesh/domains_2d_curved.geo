// Gmsh project created on Wed May 11 18:17:21 2022
SetFactory("OpenCASCADE");
//+
h = 0.05;
Mesh.CharacteristicLengthMax = h;
//+
Disk(1) = {0, 0, 0, 2, 2};
Disk(2) = {0,0.5,0,0.25,0.15};
Disk(3) = {0,-0.5,0,0.35,0.25};
b() = BooleanFragments{ Surface{1}; Delete; }{Surface{2,3}; Delete;};
Physical Surface("MatTwo") = {2,3};
Physical Surface("MatOne") = {3};
Physical Curve("RequiredBoundaryOfRequiredElements") = {1,2};
Physical Curve("OtherRequiredBoundary") = {3};