SetFactory("OpenCASCADE");
//+
h = 0.05;

Nv = DefineNumber[ {{ Nv }}, Name "Parameters/Nv" ];
Nh = DefineNumber[ {{ Nh }}, Name "Parameters/Nh" ];
L = DefineNumber[ {{ L }}, Name "Parameters/L" ];
h = DefineNumber[ {{ h }}, Name "Parameters/h" ];
hmax=0.5;
Mesh.CharacteristicLengthMax = hmax;
//+
// semi-boxes 
For r In {1:Nv}
    For s In {1:Nh}
        k = s + Nh * (r-1);
        Printf("%g, %g, box = %g", s, r, k);
        {{ ElementShape }}(k) = {{ ElementArgs }};
    EndFor
EndFor

S[]=BooleanFragments{ {{ eltDim }}{1}; Delete; }{ {{ eltDim }}{2:Nh*Nv}; Delete;};
Characteristic Length{ PointsOf{ Surface{ : }; } } = h;

//+
// Physical {{ eltDim }} (Sprintf("fin-%g",r)) = r+1;
For r In {1:Nv}
    For s In {1:Nh}
        k = s + Nh * (r-1);
        Printf("Physical Surface = %g ", k);
        Physical {{ eltDim }}(Sprintf("mat_%g",k)) = {k};
    EndFor 
EndFor 

bdy[] = CombinedBoundary { {{ eltDim }}{:}; };

Physical {{ eltDimM1 }}("Tfourier") = {4};

Physical {{ eltDimM1 }}("Tflux") =  {{ diffVal }};

For ii In { 1 : (#bdy[]-1) }
    Printf("boundary number %g = %g", ii, Abs(bdy[ii]));
    If (Abs(bdy[ii]) != {{ diffVal }})
        // Printf("boundary out number %g = %g", ii, Abs(bdy[ii]));
        Physical {{ eltDimM1 }}("Tfourier") += Abs(bdy[ii]);   
    Else
        // Printf("boundary In number %g = %g", ii, Abs(bdy[ii]));
        Physical {{ eltDimM1 }}("Tflux") = Abs(bdy[ii]);
    EndIf
EndFor

// Mesh 2;
// Show "*";
