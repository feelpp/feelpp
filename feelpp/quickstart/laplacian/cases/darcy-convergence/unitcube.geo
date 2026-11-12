SetFactory("OpenCASCADE");

h = 0.5;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

Box(1) = {0, 0, 0, 1, 1, 1};

Physical Surface("Dirichlet") = {1, 2, 3, 4, 5, 6};
Physical Volume("Omega") = {1};
