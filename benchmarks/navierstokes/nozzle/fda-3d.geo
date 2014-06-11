c=1;
hmax =0.0018;
hmin = 0.00017;
lo=0.18;

Mesh.CharacteristicLengthMin=hmin;
Mesh.CharacteristicLengthMax=hmax;
Mesh.CharacteristicLengthExtendFromBoundary=1;
Mesh.CharacteristicLengthFromPoints=1;

Point(5) =  {0.006,   0,  lo, hmax};
Point(6) =  {0,       0,  lo, hmax};
Point(7) =  {0.006,   0,  0, hmax};
Point(8) =  {0,       0,  0, hmin};
Point(9) =  {0.002,   0,  0, hmin};
Point(10) = {0,       0,  0, hmin};
Point(11) = {0.002,   0,  -0.04, hmin};
Point(12) = {0,       0,  -0.04, hmin};
Point(13) = {0.006,   0,  -0.062885, hmax};
Point(14) = {0,       0,  -0.062885, hmax};
Point(15) = {0.006,   0,  -0.182885, hmax};
Point(16) = {0,       0,  -0.182885, hmax};

Point(17) =  {0.006,   0,  -0.088, hmax};
Point(18) =  {0,   0,  -0.088, hmax};

Point(19) =  {0.006,   0,  -0.064, hmax};
Point(20) =  {0,   0,  -0.064, hmax};


Point(23) =  {0.002,   0,  -0.02, hmax};
Point(24) =  {0,   0,  -0.02, hmax};

Point(25) =  {0.002,   0,  -0.008, hmax};
Point(26) =  {0,   0,  -0.008, hmax};

Point(27) =  {0.006,   0,  0.008, hmax};
Point(28) =  {0,   0,  0.008, hmax};

Point(29) =  {0.006,   0,  0.016, hmax};
Point(30) =  {0,   0,  0.016, hmax};

Point(31) =  {0.006,   0,  0.024, hmax};
Point(32) =  {0,   0,  0.024, hmax};

Point(33) =  {0.006,   0,  0.032, hmax};
Point(34) =  {0,   0,  0.032, hmax};

Point(35) =  {0.006,   0,  0.06, hmax};
Point(36) =  {0,   0,  0.06, hmax};

Point(37) =  {0.006,   0,  0.08, hmax};
Point(38) =  {0,   0,  0.08, hmax};

Line(1) = {5, 6};
Line(2) = {8, 6};
Line(3) = {8, 10};
Line(4) = {9, 7};
Line(5) = {7, 5};
Line(6) = {9, 11};
Line(7) = {10, 12};
Line(8) = {11, 13};
Line(9) = {12, 14};
Line(10) = {14, 16};
Line(11) = {13, 15};
Line(12) = {15, 16};
Line(19) = {11, 12};
Line(20) = {9, 10};

Line(21) = {17, 18};
Line(22) = {19, 20};
Line(24) = {23, 24};
Line(25) = {25, 26};
Line(26) = {27, 28};
Line(27) = {29, 30};
Line(28) = {31, 32};
Line(29) = {33, 34};
Line(30) = {35, 36};
Line(31) = {37, 38};

Line Loop(17) = {12, -10, -9, -7, -3, 2, -1, -5, -4, 6, 8, 11};
Plane Surface(17) = {17}; 


Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {  Surface{17};}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {  Surface{76};}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {  Surface{121};}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {  Surface{166};}


Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{21};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{211};
}


Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{22};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{217};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{24};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{220};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{25};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{214};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{26};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{223};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{27};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{226};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{28};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{229};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{29};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{232};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{30};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{235};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{31};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{238};
}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{265};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{244};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{250};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{259};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{241};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{247};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{253};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{268};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{274};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{262};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{280};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{283};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{288};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{286};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{277};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{271};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{295};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{298};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{301};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{304};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{307};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{310};
}

Physical Surface("inlet") = {47,92, 137, 182};
Physical Surface("outlet") = {190, 55, 100, 145};
Physical Surface("wall") = {165, 210, 75, 120, 161, 206, 71, 116, 157, 202, 67, 112, 157, 202, 67, 112, 198, 63, 108, 153, 194, 59, 104, 149};

Physical Surface("face1") = {213, 216, 258, 234};
Physical Surface("face2") = {219, 222, 228, 246};
Physical Surface("face3") = {225, 270, 240, 291};
Physical Surface("face4") = {231, 276, 252, 294};
Physical Surface("face5") = {237, 315, 264, 297};
Physical Surface("face6") = {318, 300, 243, 282};
Physical Surface("face7") = {303, 321, 249, 285};
Physical Surface("face8") = {306, 324, 255, 288};
Physical Surface("face9") = {309, 327, 261, 279};
Physical Volume(361) = {2, 1, 4, 3};