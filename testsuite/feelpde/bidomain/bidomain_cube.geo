SetFactory("OpenCASCADE");  
h=0.1;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Box(1) = {0, 0, 0, 1, 1, 1};
Box(2) = {0.2, 0.2, 0.2, 0.4, 0.4, 0.4};
BooleanFragments{ Volume{1,2}; Delete; }{}
Physical Volume("Omega1") = {2};
Physical Volume("Omega2") = {3};
