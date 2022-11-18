SetFactory("OpenCASCADE");
//+
h = 0.05;

N = DefineNumber[ 3, Name "Parameters/N" ];
L = DefineNumber[ 1, Name "Parameters/L" ];
t = DefineNumber[ 1, Name "Parameters/t" ];
d = DefineNumber[ 0, Name "Parameters/d" ];
hmax=0.5;
Mesh.CharacteristicLengthMax = hmax;
//+
// Left boxes 
For r In {1:N}
    Printf("boxe = %g",r);
    Rectangle(r) = {-L, r*(d+t), 0, L, t, 0};
EndFor
// Middle boxes 
For r In {1:N}
    Printf("boxe = %g",r+N);
    Rectangle(r+N) = {0, r*(d+t), 0, L, t, 0};
EndFor
// Right boxes 
For r In {1:N}
    Printf("boxe = %g",r);
    Rectangle(r+2*N) = {L, r*(d+t), 0, L, t, 0};
EndFor

S[]=BooleanFragments{ Surface{1}; Delete; }{ Surface{2:3*N}; Delete;};
Characteristic Length{ PointsOf{ Surface{ : }; } } = h;

//+
// Physical Surface (Sprintf("fin-%g",r)) = r+1;
For r In {1:N}
    Printf("Physical Surface = %g %g %g ", r, r + N, r + 2*N);
    Physical Surface(Sprintf("mat_%g",r)) = {r};
    Physical Surface(Sprintf("mat_%g",r+N)) = {r+N};
    Physical Surface(Sprintf("mat_%g",r+2*N)) = {r+2*N};
EndFor 

bdy[] = CombinedBoundary { Surface{:}; };

Physical Curve("Tfourier") = {4};

Physical Curve("Tflux") =  1;

For ii In { 1 : (#bdy[]-1) }
    If (Abs(bdy[ii]) != 1)
        Printf("boundary out number %g = %g", ii, Abs(bdy[ii]));
        Physical Curve("Tfourier") += Abs(bdy[ii]);   
    Else
        Printf("boundary In number %g = %g", ii, Abs(bdy[ii]));
        Physical Curve("Tflux") = Abs(bdy[ii]);
    EndIf
EndFor

// Mesh 2;
// Show "*";
