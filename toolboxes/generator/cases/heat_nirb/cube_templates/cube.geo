SetFactory("OpenCASCADE");
//+
h = 0.05;

Nv = DefineNumber[ {{ Nv }}, Name "Parameters/Nv" ];
Nh = DefineNumber[ {{ Nh }}, Name "Parameters/Nh" ];
L = DefineNumber[ {{ L }}, Name "Parameters/L" ];
h = DefineNumber[ {{ h }}, Name "Parameters/h" ];
{% if dim == '3' -%}
d = DefineNumber[ {{ d }}, Name "Parameters/d" ];
Nd = DefineNumber[ {{ Nd }}, Name "Parameters/Nd" ];
{% endif %}
hmax = 0.1;
Mesh.CharacteristicLengthMax = hmax;

// semi-boxes
{% if dim == '2' -%}
For r In {1:Nv}
    For s In {1:Nh}
        k = s + Nh * (r-1);
        Printf("%g, %g, box = %g", s, r, k);
        {{ ElementShape }}(k) = {{ ElementArgs }};
    EndFor
EndFor
{% else %}
For z In {1:Nd}
    For y In {1:Nv}
        For x In {1:Nh}
            k = x + Nh *( (y-1) + Nv * (z-1));
            Printf("%g, %g, %g, box = %g", x, y, z, k);
            {{ ElementShape }}(k) = {{ ElementArgs }};
        EndFor
    EndFor
EndFor
{% endif %}

S[]=BooleanFragments{ {{ eltDim }}{1}; Delete; }{ {{ eltDim }}{2:Nh*Nv}; Delete;};
Characteristic Length{ PointsOf{ Surface{ : }; } } = h;

// Physical {{ eltDim }} (Sprintf("fin-%g",r)) = r+1;
{% if dim == '2' -%}
For r In {1:Nv}
    For s In {1:Nh}
        k = s + Nh * (r-1);
        Printf("Physical Surface = %g ", k);
        Physical {{ eltDim }}(Sprintf("mat_%g",k)) = {k};
    EndFor
EndFor
{% else %}
For z In {1:Nd}
    For y In {1:Nv}
        For x In {1:Nh}
            k = x + Nh *( (y-1) + Nv * (z-1));
            Printf("Physical Surface = %g ", k);
            Physical {{ eltDim }}(Sprintf("mat_%g",k)) = {k};
        EndFor
    EndFor
EndFor
{% endif %}

bdy[] = CombinedBoundary { {{ eltDim }}{:}; };

Physical {{ eltDimM1 }}("Tfourier") = {{ fourierVal }};
{% if dim == '3' -%}
Physical {{ eltDimM1 }}("Tflux") =  {};
{% else %}
Physical {{ eltDimM1 }}("Tflux") =  {{ diffVal }};
{% endif %}

For ii In { 1 : (#bdy[]-1) }
    If (Abs(bdy[ii]) != {{ diffVal }})
        Printf("boundary number %g = %g : boundary OUT", ii, Abs(bdy[ii]));
        Physical {{ eltDimM1 }}("Tfourier") += Abs(bdy[ii]);
    Else
        Printf("boundary number %g = %g : boundary IN", ii, Abs(bdy[ii]));
        Physical {{ eltDimM1 }}("Tflux") += Abs(bdy[ii]);
    EndIf
EndFor

// Mesh 2;
// Show "*";
