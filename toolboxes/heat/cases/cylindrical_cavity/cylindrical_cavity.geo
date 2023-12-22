// Gmsh project created on Fri Dec 16 10:05:59 2022

h=0.1;

Dext=0.6;
Dint=0.5;
L=1;
d=0.25; 

Z = 0.1;

SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 1, Dext, 2*Pi};
Cylinder(2) = {0, 0, 0, 0, 0, 1, Dint, 2*Pi};
//+
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Cylinder(4) = {0,0,-Z, 0, 0, Z, Dext, 2*Pi};
//+
Cylinder(5) = {0,0,-Z, 0, 0, Z, Dint, 2*Pi};
//+
BooleanFragments{ Volume{4}; Delete; }{ Volume{5}; Delete; }
//+
Cylinder(10) = {0,0,L, 0, 0, Z, Dext, 2*Pi};
//+
Cylinder(11) = {0,0,L, 0, 0, Z, d, 2*Pi};
//+
Cylinder(12) = {0,0,L, 0, 0, Z, Dint, 2*Pi};
//+
BooleanFragments{ Volume{12}; Delete; }{ Volume{11}; Delete; }//+
Dilate {{0, 0, 0}, {1, 1, L}} {
  Volume{2}; Volume{3}; 
}
//+
BooleanFragments{ Volume{10}; Delete; }{ Volume{12}; Delete; }
//+
Coherence;

//+
Cylinder(40) = {0,0,-Z*0.01, 0, 0, Z*0.01, Dint, 2*Pi};
//+
Cylinder(41) = {0,0,-Z*0.01, 0, 0, Z*0.01, Dint-0.01*Dint, 2*Pi};
//+
BooleanFragments{ Volume{40}; Delete; }{ Volume{41}; Delete; }
//+
BooleanFragments{ Volume{5}; Delete; }{ Volume{42}; Delete; }
//+
Physical Volume("InsulationBottom", 63) = {42};
//+
Coherence;

//+
Physical Volume("MaterialBottom", 69) = {41, 43};

//+
Physical Surface("CavityBottom", 71) = {43};
//+
Physical Surface("ExternalBoundaryBottom", 72) = {47};

//+
Cylinder(400) = {0,0,L, 0, 0, Z*0.01, Dint, 2*Pi};
//+
Cylinder(401) = {0,0,L, 0, 0, Z*0.01, Dint-0.01*Dint, 2*Pi};
//+
BooleanFragments{ Volume{400}; Delete; }{ Volume{401}; Delete; }

//+
Coherence;
//+
Physical Volume("MaterialTop", 407) = {405, 406};
//+
Physical Volume("InsulationTop", 408) = {402};
//+
Physical Surface("CavityTop", 409) = {64};
//+
Physical Surface("CavitySide", 410) = {65};

//+
Cylinder(500) = {0,0,L-Z*0.01, 0, 0, Z*0.01, Dint, 2*Pi};
//+
Cylinder(501) = {0,0,L-Z*0.01, 0, 0, Z*0.01, Dint+0.01*Dint, 2*Pi};
//+
BooleanFragments{ Volume{501}; Delete; }{ Volume{500}; Delete; }

//+
Cylinder(600) = {0,0,0, 0, 0, Z*0.01, Dint, 2*Pi};
//+
Cylinder(601) = {0,0,0, 0, 0, Z*0.01, Dint+0.01*Dint, 2*Pi};
//+
BooleanFragments{ Volume{601}; Delete; }{ Volume{600}; Delete; }
//+
Coherence;
//+
Physical Volume("InsulationL1", 604) = {501};
//+
Physical Volume("InsulationL2", 605) = {601};
//+
Physical Volume("LateralVolume", 606) = {603};
//+
Physical Surface("ExternalBoundaryLateral", 606) = {68};
//+
Physical Surface("ExternalBoundaryTop", 607) = {61};



Characteristic Length {PointsOf{Surface{:};}} = 0.5*h;
