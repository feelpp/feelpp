SetFactory("OpenCASCADE");
//+
h = 0.05;
P = DefineNumber[ 1, Name "Parameters/P" ];
N = DefineNumber[ 4, Name "Parameters/N" ];
L = DefineNumber[ 2.5, Name "Parameters/L" ];
R = DefineNumber[ 0.7, Name "Parameters/R" ];
t = DefineNumber[ 0.25, Name "Parameters/t" ];
d = DefineNumber[ 0.5, Name "Parameters/d" ];
hmax=0.5;
Mesh.CharacteristicLengthMax = hmax;
//+
Box(1) = {0, 0, 0, 1, P, N*(d+t)+t};
For r In {1:N}
    Printf("fin = %g",r);
    Box(r+1) = {-L, 0, r*(d+t), 2*L+1, P, t};
EndFor
S[]=BooleanFragments{ Volume{1}; Delete; }{Volume{2:N+1}; Delete;};
Characteristic Length{ PointsOf{ Surface{ : }; } } = h;
//+
// Physical Volume (Sprintf("fin-%g",r)) = r+1;
Physical Volume("Post") = {1:N*2};
i = 1;
For r In {2*N+1:#S[]:2}
    Physical Volume(Sprintf("Fin_%g",i)) = {r,r+1};
    i = i+1;
EndFor
bdy[] = CombinedBoundary { Volume{:}; };
Physical Surface("Gamma_ext") = bdy[0];
For ii In { 1 : (#bdy[]-1) }
    Printf("boundary number %g = %g", ii, bdy[ii]);
    If (Abs(bdy[ii]) != 5)
        Physical Surface("Gamma_ext") += Abs(bdy[ii]);
    Else
        Physical Surface("Gamma_root") = Abs(bdy[ii]);
    EndIf
EndFor
//Mesh 2; //+~
