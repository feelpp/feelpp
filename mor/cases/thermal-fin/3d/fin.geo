SetFactory("OpenCASCADE");
//+
h = 0.03125;

N = DefineNumber[ 4, Name "Parameters/N" ];
L = DefineNumber[ 2.5, Name "Parameters/L" ];
t = DefineNumber[ 0.25, Name "Parameters/t" ];
d = DefineNumber[ 0.75, Name "Parameters/d" ];
hmax=0.5;
Mesh.CharacteristicLengthMax = hmax;
//+
Box(1) = {0, 0, 0, 1, 1, N*(d+t)+t};
For r In {1:N}
    Printf("fin = %g",r);
    Box(r+1) = {-L, -L, r*(d+t), 2*L+1, 2*L+1, t};
EndFor

S[]=BooleanFragments{ Volume{1}; Delete; }{ Volume{2:N+1}; Delete;};
Characteristic Length{ PointsOf{ Surface{ : }; } } = h;

//+
// Physical Volume (Sprintf("fin-%g",r)) = r+1;

Physical Volume("Post") = {1:2*N};
finid = 1;
For r In {2*N+1:#S[]:1 }
    Printf("Fin number %g %g ", r, finid);
    Physical Volume(Sprintf("Fin_%g",finid)) = { r };
    finid += 1;
EndFor

bdy[] = CombinedBoundary { Volume{:}; };

Physical Surface("Gamma_ext") = bdy[0];

For ii In { 1 : (#bdy[]-1) }
    If (Abs(bdy[ii]) != 5)
        Printf("boundary number %g = %g", ii, Abs(bdy[ii]));
        Physical Surface("Gamma_ext") += Abs(bdy[ii]);   
    Else
        Physical Surface("Gamma_root") = Abs(bdy[ii]);
    EndIf
EndFor

//Mesh 2;
//Save "savedmesh.msh";
// Show "*";
