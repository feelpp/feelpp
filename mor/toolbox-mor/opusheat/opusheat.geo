m=1;mm=10^-3;

//
// Integrated circuit : IC
//
// thickness
e_IC  = 2*2e-3*m;
// length
L_IC  = 2e-2*m;
// position of the first IC
h_1   = 20*mm;
// position of the second IC
h_2   = 70*mm;


//
// PCB
//
// thickness
e_PCB = 2*2e-3*m;
// height
h_PCB = 13e-2*m;


//
// Air
//
// thickness
e_A = 5e-2*m;
h=0.2*mm - 1e-8;
h=0.001;
hs=h;
p1=newp;Point(p1) = {0,0,0,hs};
p2=newp;Point(p2) = {e_PCB,0,0,hs};

p3=newp;Point(p3) = {e_PCB,h_1,0,hs};
p4=newp;Point(p4) = {e_PCB,h_1+L_IC,0,hs};

p5=newp;Point(p5) = {e_PCB,h_2,0,h};
p6=newp;Point(p6) = {e_PCB,h_2+L_IC,0,hs};

p7=newp;Point(p7) = {e_PCB,h_PCB,0,hs};
p8=newp;Point(p8) = {0,h_PCB,0,hs};

air_p1=p2;
air_p21=newp;Point(air_p21) = {e_PCB+e_IC,0,0,hs};
air_p22=newp;Point(air_p22) = {e_PCB+e_A,0,0,hs};
air_p31=newp;Point(air_p31) = {e_PCB+e_A,h_PCB,0,hs};
air_p32=newp;Point(air_p32) = {e_PCB+e_IC,h_PCB,0,hs};
air_p4=p7;

air_p5=p3;
air_p51=newp;Point(air_p51) = {e_PCB+e_IC,h_1,0,hs};
air_p6=p4;
air_p61=newp;Point(air_p61) = {e_PCB+e_IC,h_1+L_IC,0,hs};

air_p7=p5;
air_p71=newp;Point(air_p71) = {e_PCB+e_IC,h_2,0,hs};
air_p8=p6;
air_p81=newp;Point(air_p81) = {e_PCB+e_IC,h_2+L_IC,0,hs};

ic1_p1=p3;
ic1_p2=air_p51;
ic1_p3=air_p61;
ic1_p4=p4;

ic2_p1=p5;
ic2_p2=air_p71;
ic2_p3=air_p81;
ic2_p4=p6;

Line(1) = {1, 2};
Line(2) = {2, 9};
Line(3) = {9, 10};
Line(4) = {10, 11};
Line(5) = {11, 12};
Line(6) = {12, 7};
Line(7) = {7, 8};
Line(8) = {7, 8};
Line(9) = {8, 1};
Line(10) = {2, 3};
Line(11) = {3, 4};
Line(12) = {4, 5};
Line(13) = {5, 6};
Line(14) = {6, 7};
Line(15) = {9, 13};
Line(16) = {13, 14};
Line(17) = {14, 15};
Line(18) = {15, 16};
Line(19) = {16, 12};
Line(20) = {3, 13};
Line(21) = {4, 14};
Line(22) = {5, 15};
Line(23) = {6, 16};
Line Loop(24) = {4, 5, -19, -18, -17, -16, -15, 3};
Plane Surface(25) = {24};
Line Loop(26) = {19, 6, -14, 23};
Plane Surface(27) = {26};
Line Loop(28) = {18, -23, -13, 22};
Plane Surface(29) = {28};
Line Loop(30) = {17, -22, -12, 21};
Plane Surface(31) = {30};
Line Loop(32) = {16, -21, -11, 20};
Plane Surface(33) = {32};
Line Loop(34) = {10, 20, -15, -2};
Plane Surface(35) = {34};
Line Loop(36) = {9, 1, 10, 11, 12, 13, 14, 7};
Plane Surface(37) = {36};

//Physical Line(38) = {10, 11, 12, 13, 14};
Physical Line("Gamma_1") = {9};
Physical Line("Gamma_2") = {4};
Physical Line("Gamma_IC1_PCB") = {11};
Physical Line("Gamma_IC2_PCB") = {13};

Physical Line("Gamma_4_AIR1") = {2};
Physical Line("Gamma_4_AIR4") = {3};
Physical Line("Gamma_3_AIR4") = {5};
Physical Line("Gamma_3_AIR3") = {6};
Physical Line("Gamma_3_PCB") = {7};
Physical Line("Gamma_4_PCB") = {1};
Physical Surface("PCB") = {37};
Physical Surface("AIR") = {35, 31, 27,25};
Physical Surface("IC1") = {33};
Physical Surface("IC2") = {29};
