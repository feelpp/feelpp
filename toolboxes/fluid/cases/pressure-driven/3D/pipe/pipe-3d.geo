SetFactory("OpenCASCADE");

height = 1;
r = 0.2;
h=0.1;
Cylinder(1) = {0, 0, 0, height, 0, 0, r, 2*Pi};
Mesh.CharacteristicLengthFromCurvature=10;
Characteristic Length{ PointsOf{ Surface{1}; } } = h;
Physical Surface("wall") = {1};
Physical Surface("outlet") = {2};
Physical Surface("inlet") = {3};

Physical Volume("Omega") = {1};