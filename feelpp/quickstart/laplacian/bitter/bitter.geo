h = 2;

part = DefineNumber[36, Name "Bitter/angle : Pi on"];
angle = Pi/part;
Ri = DefineNumber[305, Name "Bitter/Ri"];
Re = DefineNumber[401, Name "Bitter/Re"];
Rc1 = DefineNumber[319.26, Name "Cooling/Rc1"];
Rc2 = DefineNumber[349.71, Name "Cooling/Rc2"];
Rc3 = DefineNumber[383.11, Name "Cooling/Rc3"];
Rt = DefineNumber[5, Name "Cooling/Rt"];

Lc = DefineNumber[5.9, Name "Cooling/Lc"];
Wc = DefineNumber[1.1, Name "Cooling/Wc"];

height = DefineNumber[1, Name "Bitter/Lz"];

// Rectangle for cooling holes
Point(1) = {0,0,0};
Point(2) = {0,Lc/2,0,h};
Point(3) = {-Wc/2,Lc/2,0,h};
Point(4) = {0,(Lc+Wc)/2,0,h};
Point(5) = {Wc/2,Lc/2,0,h};
Point(6) = {0,-Lc/2,0,h};
Point(7) = {-Wc/2,-Lc/2,0,h};
Point(8) = {0,-(Lc+Wc)/2,0,h};
Point(9) = {Wc/2,-Lc/2,0,h};
Point(10) = {Wc/2,0,0,h};
Point(11) = {-Wc/2,0,0,h};

Circle(1) = {3,2,4};
Circle(2) = {4,2,5};
Line(3) = {5,10};
Line(4) = {10,9};
Circle(5) = {9,6,8};
Circle(6) = {8,6,7};
Line(7) = {7,11};
Line(8) = {11,3};

Rc[] = {Rc1,Rc3};
a[] = {1,3,-1,-3};
For i In {0:#Rc[]-1}
    For j In {0:#a[]-1}
        out[] = Rotate {{0,0,1}, {0,0,0}, a[j]*angle/4} {Translate {Rc[i],0,0} {Duplicata {Line{1,2,3,4,5,6,7,8}; } } };
        cLL[] += newll; Line Loop(newll) = {out[]};
    EndFor
EndFor
a[] = {-1,1};
For j In {0:#a[]-1}
    out[] = Rotate {{0,0,1}, {0,0,0}, a[j]*angle/2} {Translate {Rc2,0,0} {Duplicata {Line{1,2,3,4,5,6,7,8}; } } };
    cLL[] += newll; Line Loop(newll) = {out[]};
EndFor

// Half Rectangle Points
bp[] = Rotate {{0,0,1}, {0,0,0}, angle} {Translate {Rc2,0,0} { Duplicata {Point{6,7,8,9,10,11}; } } };
bp[] +=Rotate {{0,0,1}, {0,0,0}, -angle} {Translate {Rc2,0,0} { Duplicata {Point{2,3,4,5,10,11}; } } };

Recursive Delete {Line{1,2,3,4,5,6,7,8}; }

// Tie Rodes
cTc = newp; Point(cTc) = {Rc2,0,0,1.2*h};
rTc1 = newp; Point(rTc1) = {Rc2,Rt,0,1.2*h};
rTc2 = newp; Point(rTc2) = {Rc2+Rt,0,0,1.2*h};
rTc3 = newp; Point(rTc3) = {Rc2,-Rt,0,1.2*h};
rTc4 = newp; Point(rTc4) = {Rc2-Rt,0,0,1.2*h};

cTc1 = newl; Circle(cTc1) = {rTc1,cTc,rTc2};
cTc2 = newl; Circle(cTc2) = {rTc2,cTc,rTc3};
cTc3 = newl; Circle(cTc3) = {rTc3,cTc,rTc4};
cTc4 = newl; Circle(cTc4) = {rTc4,cTc,rTc1};

ll = newll; Line Loop(newll) = {cTc1,cTc2,cTc3,cTc4};

// Bitter
R[] = {Ri,Re};
a[] = {-1,1};
For i In {0:#R[]-1}
    For j In {0:#a[]-1}
        p[] += newp; Point(newp) = {R[i]*Cos(a[j]*angle),R[i]*Sin(a[j]*angle),0,1.5*h};
    EndFor
EndFor

bL[]  = newl; Circle(newl) = {p[0],1,p[1]};
bL[] += newl; Line(newl) = {p[1],bp[5]};
bL[] += newl; Line(newl) = {bp[5],bp[1]};
bL[] += newl; Circle(newl) = {bp[1],bp[0],bp[2]};
bL[] += newl; Circle(newl) = {bp[2],bp[0],bp[3]};
bL[] += newl; Line(newl) = {bp[3],bp[4]};
bL[] += newl; Line(newl) = {bp[4],p[3]};
bL[] += newl; Circle(newl) = {p[3],1,p[2]};
bL[] += newl; Line(newl) = {p[2],bp[10]};
bL[] += newl; Line(newl) = {bp[10],bp[9]};
bL[] += newl; Circle(newl) = {bp[9],bp[6],bp[8]};
bL[] += newl; Circle(newl) = {bp[8],bp[6],bp[7]};
bL[] += newl; Line(newl) = {bp[7],bp[11]};
bL[] += newl; Line(newl) = {bp[11],p[0]};

bLL = newll; Line Loop(bLL) = {bL[]};

bS = news; Plane Surface(bS) = {bLL, ll, cLL[]};

bV[] = Extrude {0,0,height } {Surface{bS}; };

cS[] = {bV[4],bV[5],bV[6],bV[7],bV[11],bV[12],bV[13],bV[14]};
For i In {20:#bV[]-1}
    cS[] += bV[i];
EndFor


Physical Volume("omega") = {bV[1]};
//Robin
Physical Surface("Robin") = {bV[2],bV[9],bV[16],bV[17],bV[18],bV[19],cS[]};
//Physical Surface("Channel0") = {bV[2]};
//Physical Surface("Channel1") = {bV[9]};
//Physical Surface("TieRods") = {bV[16],bV[17],bV[18],bV[19]};
//Physical Surface("CoolingHoles") = {cS[]};

// Neumann
Physical Surface("Neumann") = {bV[0],bS,bV[3],bV[8],bV[10],bV[15]};
//Physical Surface("top") = {bV[0]};
//Physical Surface("bottom") = {bS};
//Physical Surface("V0") = {bV[3],bV[8]};
//Physical Surface("V1") = {bV[10],bV[15]};
