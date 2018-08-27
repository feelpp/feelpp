m=1;mm=10^-3;

//
// Integrated circuit : IC
//
// thickness
e_IC  = 2e-3*m;
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
e_PCB = 2e-3*m;
// height
h_PCB = 13e-2*m;


//
// Air
//
// thickness
e_A = 5e-2*m;
h=0.2*mm - 1e-8;

Include "geometry_heat.geo";
