h=0.3;

hasAir = 1;
hasAir = DefineNumber[1, Name "Air/enable", Choices{0="false", 1="true"}];

r1=1;
r1 = DefineNumber[1, Name "conductor/internal radius", Min 0.1];
r2=2;
r2 = DefineNumber[2, Name "conductor/external radius", Min r1+0.1];
L=r2;
L = DefineNumber[r2, Name "conductor/height", Min 0.1];

factor = 0.5;
factor = DefineNumber[0.5, Name "conductor/angle (*Pi)", Min 0.01, Max 1.99];
angle = factor*Pi;

rext=4*r2;
rext = DefineNumber[4*r2, Name "Air/external radius", Min r2+0.1, Visible hasAir];
rinf=6*r2;
rinf = DefineNumber[6*r2, Name "Air/infinite radius", Min rext+0.1, Visible hasAir];
hext=5*h;
hinf=10*h;

omega[] = {};
bottom[] = {};
Rext[] = {};
top[] = {};
Rint[] = {};
air[] = {};
border[] = {};
plan[] = {};

// Conductor
Point(0) = {0,0,0,h};
Point(1) = {r1,0,-L/2,h};
Point(2) = {r2,0,-L/2,h};
Point(3) = {r2,0,L/2,h};
Point(4) = {r1,0,L/2,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

If ( hasAir == 1 )
    Point(5) = {0,0,-rext,hext};
    Point(6) = {rext,0,0,hext};
    Point(7) = {0,0,rext,hext};

    Circle(5) = {5,0,6};
    Circle(6) = {6,0,7};
    Line(7) = {7,5};

    Line Loop(2) = {5,6,7};
    Plane Surface(2) = {2,1};

    Point(8) = {0,0,-rinf,hinf};
    Point(9) = {rinf,0,0,hinf};
    Point(10) = {0,0,rinf,hinf};

    Circle(8) = {8,0,9};
    Line(9) = {9,6};
    Line(10) = {5,8};

    Line Loop(3) = {8,9,-5,10};
    Plane Surface(3) = {3};
    
    Circle(11) = {9,0,10};
    Line(12) = {10,7};

    Line Loop(4) = {11,12,-6,-9};
    Plane Surface(4) = {4};

    plan[] = {2,3,4};
    Physical Surface("OX") = {plan[]};
EndIf

outS = 1;
Physical Surface("out") = {1};

For i In {0:3}
    out[] = Extrude{{0,0,1}, {0,0,0}, angle/4} {Surface{outS,plan[]};};
    outS = out[0];
    omega[] += {out[1]};
    bottom[] += {out[2]};
    Rext[] += {out[3]};
    top[] += {out[4]};
    Rint[] += {out[5]};
    If ( hasAir == 1 )
        air[] += {out[7],out[15],out[20]};
        border[] += {out[16],out[21]};
        plan[] = {out[6],out[14],out[19]};
    EndIf
EndFor

Physical Volume("omega") = {omega[]};
Physical Surface("in") = {outS};
Physical Surface("bottom") = {bottom[]};
Physical Surface("Rext") = {Rext[]};
Physical Surface("top") = {top[]};
Physical Surface("Rint") = {Rint[]};

If ( hasAir == 1 )
    Physical Volume("air") = {air[]};
    Physical Surface("OY") = {plan[]};
    Physical Surface("Border") = {border[]};
EndIf
