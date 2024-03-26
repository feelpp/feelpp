// Gmsh project created on Fri Mar 18 16:22:09 2016

SetFactory("OpenCASCADE");  
h=0.1;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
// using opencascade kernel create a rectangle inside another one 
Rectangle(1) = {0, 0, 0, 1, 1, 0};
Rectangle(2) = {0.2, 0.2, 0, 0.4, 0.4, 0};
BooleanFragments{ Surface{1,2}; Delete; }{}
//Characteristic Length{ PointsOf{ Surface{:}; } } = h;
Physical Line("Gamma1_D") = {1,2};
Physical Line("Gamma1_N") = {3,4};
Physical Line("Gamma2_D") = {5,6};
Physical Line("Gamma2_N") = {7,8};
Physical Surface("Omega1") = {2};
Physical Surface("Omega2") = {3};
