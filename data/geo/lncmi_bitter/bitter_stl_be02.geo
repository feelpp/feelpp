Mesh.RemeshParametrization=1; //(0) harmonic (1) conformal
Mesh.RemeshAlgorithm=1; //(0) nosplit (1) automatic (2) split metis

Mesh.Algorithm=6;
Mesh.Algorithm3D=4;
//Mesh.Partitioner=1;
//Mesh.NbPartitions=32;

//Mesh.CharacteristicLengthFactor=0.05;
//Mesh.CharacteristicLengthFactor=0.2;

//Mesh.CharacteristicLengthMax = 1.0;

Merge "BE-02-01M.stp";

Physical Volume("omega") = {1};

Physical Surface("V0") = {1, 53};
Physical Surface("V1") = {56, 57};
Physical Surface("TieRods") = {51, 52};
Physical Surface("Channel0") = {54, 55};
Physical Surface("Channel1") = {7, 58};
Physical Surface("CoolingHoles") = {3, 6, 5, 35, 38, 37, 36, 32, 33, 8, 34, 31, 10, 9, 18, 17, 30, 29, 22, 21, 26, 25, 27, 28, 15, 16, 19, 20, 23, 24, 50, 47, 48, 14, 11, 12, 39, 40, 43, 44, 46, 42, 49, 13, 41, 45};
Physical Surface("top") = {2};
Physical Surface("bottom") = {4};
