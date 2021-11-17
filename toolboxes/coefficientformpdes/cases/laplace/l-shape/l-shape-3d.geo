SetFactory("OpenCASCADE");

h = 0.1; // or whatever
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

Box(1)= {0,0,0,1,1,1};
Box(2)= {0.5,0.5,0.5,0.5,0.5,0.5};
V[]=BooleanDifference{Volume{1}; Delete; }{Volume{2}; Delete;} ;
b[]=Boundary{Volume{V[]};};
Physical Surface("Gamma")={b[]};
Physical Volume("Omega")={V[]};
