SetFactory("OpenCASCADE");
//+
hs = 0.005;

Rectangle(1) = {-0.2, -0.2, 0, 1.4, 1.4, 0};
Rectangle(2) = {0, 0, 0, 1, 1, 0};

S[]=BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete;};//+
Physical Curve("Gamma_1") = {8};
//+
Physical Curve("Gamma_2") = {7};
//+
Physical Curve("Gamma_3") = {6};
//+
Physical Curve("Gamma_4") = {5};

Physical Surface("Surface") = {1};

Characteristic Length { PointsOf{ Surface{1}; }} = hs;