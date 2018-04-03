h = 1e-07;

Lx = DefineNumber[480*10^(-9), Name "x Lenght"];
Ly = DefineNumber[320*10^(-9), Name "y Lenght"];
Lz = DefineNumber[320*10^(-9), Name "z Lenght"];
Lj = DefineNumber[160*10^(-9), Name "j Lenght"];
t0x = DefineNumber[10*10^(-9), Name "Gate Lenght"];

// Discretization
npt = DefineNumber[20, Name "Nb point interp"];
l = DefineNumber[3.5, Name "Interv for tanh"];

A = newp; Point(A) = {0,0,0,h};
B = newp; Point(B) = {Lx,0,0,h};
C = newp; Point(C) = {Lx,Ly,0,h};
D = newp; Point(D) = {2./3.*Lx,Ly,0,h/3};
E = newp; Point(E) = {2./3.*Lx-Lx/8.,Ly,0,h/3};
F = newp; Point(F) = {1./3.*Lx+Lx/8.,Ly,0,h/3};
G = newp; Point(G) = {1./3.*Lx,Ly,0,h/3};
H = newp; Point(H) = {0,Ly,0,h};
K = newp; Point(K) = {2./3.*Lx-Lx/8.,Ly+t0x,0,h/3};
L = newp; Point(L) = {1./3.*Lx+Lx/8.,Ly+t0x,0,h/3};
M = newp; Point(M) = {0,Ly-Lj,0,h};
N = newp; Point(N) = {1./3.*Lx,Ly-Lj,0,h};
P = newp; Point(P) = {Lx,Ly-Lj,0,h};
Q = newp; Point(Q) = {2./3.*Lx,Ly-Lj,0,h};

AB = newl; Line(AB) = {A,B};
BP = newl; Line(BP) = {B,P};
PC = newl; Line(PC) = {P,C};
CD = newl; Line(CD) = {C,D};
DE = newl; Line(DE) = {D,E};
EF = newl; Line(EF) = {E,F};
FG = newl; Line(FG) = {F,G};
KL = newl; Line(KL) = {K,L};
GH = newl; Line(GH) = {G,H};
HM = newl; Line(HM) = {H,M};
MA = newl; Line(MA) = {M,A};
MN = newl; Line(MN) = {M,N};
PQ = newl; Line(PQ) = {P,Q};

pListNF[0] = N;
pListQE[0] = Q;
pListGL[0] = G;
pListDK[0] = D;
For i In {1:npt-1}
    xhi = i*Lx/8./npt;
    eta = 64*Lj/(Lx*Lx)*xhi*xhi;
    xNF = 1./3.*Lx + xhi;
    xQE = 2./3.*Lx - xhi;
    y = eta + Ly - Lj;
    pListNF[i] = newp;
    Point(pListNF[i]) = {xNF,y,0,h};
    pListQE[i] = newp;
    Point(pListQE[i]) = {xQE,y,0,h};
    eta = t0x/2.*(Tanh(16.*l/Lx*xhi-l)+1);
    xGL = xNF;
    xDK = xQE;
    y = Ly + eta;
    pListGL[i] = newp;
    Point(pListGL[i]) = {xGL,y,0,h/3};
    pListDK[i] = newp;
    Point(pListDK[i]) = {xDK,y,0,h/3};
EndFor
pListNF[npt] = F;
pListQE[npt] = E;
pListGL[npt] = L;
pListDK[npt] = K;

NF = newl; Spline(NF) = pListNF[];
QE = newl; Spline(QE) = pListQE[];
GL = newl; Spline(GL) = pListGL[];
DK = newl; Spline(DK) = pListDK[];

Bll = newll; Line Loop(Bll) = {AB,BP,PQ,QE,EF,-NF,-MN,MA};
OB = news; Plane Surface(OB) = {Bll};

Dll = newll; Line Loop(Dll) = {PC,CD,DE,-QE,-PQ};
OD = news; Plane Surface(OD) = {Dll};

Sll = newll; Line Loop(Sll) = {FG,GH,HM,MN,NF};
OS = news; Plane Surface(OS) = {Sll};

Gll = newll; Line Loop(Gll) = {DE,EF,FG,GL,-KL,-DK};
OG = news; Plane Surface(OG) = {Gll};

// 3D
out[] = Extrude{0,0,Lz} {Surface{OB,OD,OS,OG};};

// For i In {0:#out[]-1}
//     Printf("out[%g] : %g",i, out[i]);
// EndFor

Physical Surface("Bulk") = {out[2]};
Physical Surface("Drain") = {out[13]};
Physical Surface("Source") = {out[20]};
Physical Surface("Gate") = {out[30]};
Physical Surface("IntBD") = {out[4],out[5]};
Physical Surface("IntBG") = {out[6]};
Physical Surface("IntBS") = {out[7],out[8]};
Physical Surface("IntDG") = {out[14]};
Physical Surface("IntSG") = {out[19]};
Physical Surface("WallB") = {OB,out[0],out[3],out[9]};
Physical Surface("WallD") = {OD,out[10],out[12]};
Physical Surface("WallS") = {OS,out[17],out[21]};
Physical Surface("WallG") = {OG,out[24],out[29],out[31]};

Physical Volume("OmegaB") = {out[1]};
Physical Volume("OmegaD") = {out[11]};
Physical Volume("OmegaS") = {out[18]};
Physical Volume("OmegaG") = {out[25]};
