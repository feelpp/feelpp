SetFactory("OpenCASCADE");
//+
h = 0.05;

P = DefineNumber[ {{ P }}, Name "Parameters/P" ];
N = DefineNumber[ {{ N }}, Name "Parameters/N" ];
L = DefineNumber[ {{ L }}, Name "Parameters/L" ];
t = DefineNumber[ {{ t }}, Name "Parameters/t" ];
d = DefineNumber[ {{ d }}, Name "Parameters/d" ];
hmax=0.5;
Mesh.CharacteristicLengthMax = hmax;
//+
{{ PostShape }}(1) = {{ PostArgs }};
For r In {1:N}
    Printf("fin = %g",r);
    {{ FinShape }}(r+1) = {{ FinArgs }};
EndFor

S[]=BooleanFragments{ {{ eltDim }}{1}; Delete; }{ {{ eltDim }}{2:N+1}; Delete;};
Characteristic Length{ PointsOf{ Surface{ : }; } } = h;

//+
// Physical {{ eltDim }} (Sprintf("fin-%g",r)) = r+1;

Physical {{ eltDim }}("Post") = {1:2*N};
finid = 1;
For r In {2*N+1:#S[]:{{ step }} }
    Printf("Fin number %g %g ", r, finid);
    Physical {{ eltDim }}(Sprintf("Fin_%g",finid)) = {{ physicalArg }};
    finid += 1;
EndFor

bdy[] = CombinedBoundary { {{ eltDim }}{:}; };

Physical {{ eltDimM1 }}("Gamma_ext") = bdy[0];

For ii In { 1 : (#bdy[]-1) }
    If (bdy[ii] != {{ diffVal }})
        Printf("boundary number %g = %g", ii, Abs(bdy[ii]));
        Physical {{ eltDimM1 }}("Gamma_ext") += Abs(bdy[ii]);   
    Else
        Physical {{ eltDimM1 }}("Gamma_root") = bdy[ii];
    EndIf
EndFor

// Mesh 2;
// Show "*";
