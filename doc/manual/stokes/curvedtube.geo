Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthExtendFromBoundary=1;
Mesh.CharacteristicLengthFromPoints=1;
Mesh.ElementOrder=1;
Mesh.SecondOrderIncomplete = 0;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 4;
//Mesh.OptimizeNetgen=1;
// partitioning data
Mesh.Partitioner=1;
Mesh.NbPartitions=1;
Mesh.MshFilePartitioned=0;

h=0.1;
Point(1) = {0, 0, 0, h};
a=0.785398163397448; // pi/4
Point(2) = {a, a, 0, h};
Point(3) = {2*a, 2*a, 0, h };
Point(4) = {-a, a, 0, h};
Point(5) = {-2*a, 2*a, 0, h};

Circle(1) = {2, 1, 4};
Circle(2) = {3, 1, 5};
Line(3) = {2, 3};
Line(4) = {4, 5};
Line Loop(5) = {3, 2, -4, -1};
Plane Surface(6) = {5};
Physical Line("wall") = {1, 2};
Physical Line("outlet") = {3};
Physical Line("inlet") = {4};
Physical Surface(7) = {6};
