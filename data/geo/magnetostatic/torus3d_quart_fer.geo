Unit = 1.e-3;
lc = 0.5*Unit;

//N_Layers = 20;

r1=61.2*0.5*Unit;
r2=106.4*0.5*Unit;
h=4.61*Unit/2.;
hf=2*h;

r_inf=200*Unit;
r_ext=125*Unit;

// Define torus section
P0= newp; Point(P0) = {0,0,0, lc};
P1= newp; Point(P1) = {r1,0,-h, lc};
P2= newp; Point(P2) = {r2,0,-h, lc};
P3= newp; Point(P3) = {r2,0, h, lc};
P4= newp; Point(P4) = {r1,0, h, lc};

// Define iron section
PF1= newp; Point(PF1) = {r1/4.0,0,-hf, lc};
PF2= newp; Point(PF2) = {r2/4.0,0,-hf, lc};
PF3= newp; Point(PF3) = {r2/4.0,0, hf, lc};
PF4= newp; Point(PF4) = {r1/4.0,0, hf, lc};

L12=newl; Line(L12) = {P1, P2};
L23=newl; Line(L23) = {P2, P3};
L34=newl; Line(L34) = {P3, P4};
L41=newl; Line(L41) = {P4, P1};

LF12=newl; Line(LF12) = {PF1, PF2};
LF23=newl; Line(LF23) = {PF2, PF3};
LF34=newl; Line(LF34) = {PF3, PF4};
LF41=newl; Line(LF41) = {PF4, PF1};

L=newl; Line Loop(L) = {L12, L23, L34, L41};
S=newreg; Plane Surface(S) = {L};

LF=newl; Line Loop(LF) = {LF12, LF23, LF34, LF41};
SF=newreg; Plane Surface(SF) = {LF};


// Define ext section
P5=newp; Point(P5) = {0,0, -r_ext, 10*lc};
P6=newp; Point(P6) = {r_ext, 0, 0, 10*lc};
P7=newp; Point(P7) = {0, 0, r_ext, 10*lc};

L05=newl; Line(L05) = {P0, P5};
C56=newl; Circle(C56) = {P5, P0, P6};
C67=newl; Circle(C67) = {P6, P0, P7};
L70=newl; Line(L70) = {P7, P0};

L_ext=newl; Line Loop(L_ext) = {L05, C56, C67, L70};
S_ext=newreg; Plane Surface(S_ext) = {L_ext, -L, -LF};

// Define inf section
P8=newp; Point(P8) = {0,     0, -r_inf, 60*lc};
P9=newp; Point(P9) = {r_inf, 0, 0, 60*lc};
P10=newp; Point(P10) = {0,   0, r_inf, 60*lc};

L58=newl; Line(L58) = {P5, P8};
C89=newl; Circle(C89) = {P8, P0, P9};
C910=newl; Circle(C910) = {P9, P0, P10};
L107=newl; Line(L107) = {P10, P7};

L_inf=newl; Line Loop(L_inf) = {L58, C89, C910, L107, -C67, -C56};
S_inf=newreg; Plane Surface(S_inf) = {L_inf};

// Build 3D geom
Extrude { {0,0,1} , {0,0,0} , Pi/2. } { 
  Surface{S, S_ext, S_inf}; 
}

Extrude { {0,0,1} , {0,0,0} , Pi/2. } { 
  Surface{SF}; 
}

//  Define Physical

Physical Volume("Torus") = {1}; // Tore
Physical Volume("Air") = {2}; // Air
Physical Volume("Inf") = {3}; // Infini
Physical Volume("Iron") = {4}; // Fer

Physical Surface("V0") = {10};   // V0
Physical Surface("V1") = {46};  // V1
Physical Surface("Rint") = {45};  // Rint
Physical Surface("Rext") = {37};  // Rext
Physical Surface("HChannel1") = {33};  // Hchannel Bas
Physical Surface("HChannel2") = {41};  // Hchannel Haut

// Fer
Physical Surface("V0_Iron") = {12};   // V0
Physical Surface("V1_Iron") = {144};  // V1
Physical Surface("Rint_Iron") = {99};  // Rint
Physical Surface("Rext_Iron") = {91};  // Rext
Physical Surface("HChannel1_Iron") = {87};  // Hchannel Bas
Physical Surface("HChannel2_Iron") = {95};  // Hchannel Haut

Physical Surface("Border0") = {18, 24};  // Sym
Physical Surface("Border1") = {100, 122};  // Sym

Physical Surface("BorderInf") = {111, 114};  // Inf

Physical Line("Axis") = {7,10,13,16}; //Axis
