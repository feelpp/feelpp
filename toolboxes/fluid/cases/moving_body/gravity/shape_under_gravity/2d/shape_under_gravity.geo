SetFactory("OpenCASCADE");

h = 0.05; // or whatever
//Mesh.CharacteristicLengthMin = h/10;
//Mesh.CharacteristicLengthMax = h;

Rectangle(1) = {-5,0,0,10,10};
Rectangle(2) = {-0.05,8,0,0.1,0.2};
Characteristic Length{ PointsOf{ Line{1:4}; } } = 2*h;
Characteristic Length{ PointsOf{ Line{5:8}; } } = h/20;

S[]=BooleanFragments{Surface{1,2}; Delete; }{} ;
b[]=Boundary{Surface{S[]};};
For r In {0:#S[]-1}
    Printf("surface: %g",S[r]);
EndFor
For r In {0:#b[]-1}
    Printf("curve: %g",b[r]);
EndFor
Physical Curve("Wall")={1:4};
Physical Curve("Wall_Body")={5:8};
Physical Surface("Fluid")={3};
Physical Surface("Body")={2};