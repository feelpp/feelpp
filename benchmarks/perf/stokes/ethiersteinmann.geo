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
h = 0.1;
Mesh.RecombinationAlgorithm=0;
a=-1;
b=1;
c=-1;
d=1;
e=-1;
f=1;
Point(1) = {a,c,e,h};
Point(2) = {b,c,e,h};
Point(3) = {b,d,e,h};
Point(4) = {a,d,e,h};
Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,1};
Line(4) = {1,2};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};

Extrude Surface {6, {0,0,f-e} } {
  Layers { {(f-e)/h}, {1.0} };
};
Physical Surface("Neumann") = {6,19,27,28};
Physical Surface("Dirichlet") = {15,23};
Physical Volume("Mat1") = {1};
