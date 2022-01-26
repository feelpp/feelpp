SetFactory("OpenCASCADE");
//+
h = 0.01; 

N = DefineNumber[ 4, Name "Parameters/N" ];
L = DefineNumber[ 2.5, Name "Parameters/L" ];
t = DefineNumber[ 0.25, Name "Parameters/t" ];
d = DefineNumber[ 0.75, Name "Parameters/d" ];
hmax=0.5;
Mesh.CharacteristicLengthMax = hmax;
//+
Rectangle(1) = {0, 0, 0, 1, N*(d+t)+t, 0};
For r In {1:N}
    Printf("fin = %g",r);
    Rectangle(r+1) = {-L, r*(d+t), 0, 2*L+1, t, 0};
EndFor

S[]=BooleanFragments{ Surface{1}; Delete; }{ Surface{2:N+1}; Delete;};
Characteristic Length{ PointsOf{ Surface{ : }; } } = h;

//+
// Physical Surface (Sprintf("fin-%g",r)) = r+1;

Physical Surface("Post") = {1:2*N};
finid = 1;
For r In {2*N+1:#S[]:2 }
    Printf("Fin number %g %g ", r, finid);
    Physical Surface(Sprintf("Fin_%g",finid)) = {r,r+1};
    finid += 1;
EndFor

bdy[] = CombinedBoundary { Surface{:}; };

Physical Curve("Gamma_ext") = bdy[0];

For ii In { 1 : (#bdy[]-1) }
    If (bdy[ii] != 4)
        Printf("boundary number %g = %g", ii, Abs(bdy[ii]));
        Physical Curve("Gamma_ext") += Abs(bdy[ii]);   
    Else
        Physical Curve("Gamma_root") = bdy[ii];
    EndIf
EndFor

// Mesh 2;
// Show "*";
