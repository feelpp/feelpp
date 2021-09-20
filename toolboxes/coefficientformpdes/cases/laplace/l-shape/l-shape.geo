SetFactory("OpenCASCADE");

h = 0.2; // or whatever
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

Rectangle(1)= {0,0,0,1,2};
Rectangle(2)= {0,0,0,2,1};
S[]=BooleanFragments{Surface{1,2}; Delete; }{} ;
b[]=Boundary{Surface{S[]};};
For r In {0:#S[]-1}
    Printf("surface: %g",S[r]);
EndFor
For r In {0:#b[]-1}
    Printf("curve: %g",b[r]);
EndFor
Physical Curve("Gamma")={1,4:10};
Physical Surface("Omega")={S[]};
