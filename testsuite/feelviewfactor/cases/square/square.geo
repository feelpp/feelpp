SetFactory("OpenCASCADE");
//+
h = 0.05;

Rectangle(1) = {-0.2, -0.2, 0, 1.4, 1.4, 0};
Rectangle(2) = {0, 0, 0, 1, 1, 0};

S[]=BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete;};//+
Physical Curve("Gamma_1", 9) = {8};
//+
Physical Curve("Gamma_2", 10) = {7};
//+
Physical Curve("Gamma_3", 11) = {6};
//+
Physical Curve("Gamma_4", 12) = {5};
//+
Show "*";