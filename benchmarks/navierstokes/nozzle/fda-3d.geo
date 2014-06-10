c=1;
hmax =0.0018;
hmin = 0.0009;
lo=0.18;

Mesh.CharacteristicLengthMin=hmin;
Mesh.CharacteristicLengthMax=hmax;
Mesh.CharacteristicLengthExtendFromBoundary=1;
Mesh.CharacteristicLengthFromPoints=1;

Point(5) =  {0.006,   0,  lo, hmax};
Point(6) =  {0,       0,  lo, hmax};
Point(7) =  {0.006,   0,  0, hmin};
Point(8) =  {0,       0,  0, hmin};
Point(9) =  {0.002,   0,  0, hmin};
Point(10) = {0,       0,  0, hmin};
Point(11) = {0.002,   0,  -0.04, hmin};
Point(12) = {0,       0,  -0.04, hmin};
Point(13) = {0.006,   0,  -0.062885, hmax};
Point(14) = {0,       0,  -0.062885, hmax};
Point(15) = {0.006,   0,  -0.182885, hmax};
Point(16) = {0,       0,  -0.182885, hmax};

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
Line Loop(17) = {12, -10, -9, -7, -3, 2, -1, -5, -4, 6, 8, 11};
Plane Surface(17) = {17};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {  Surface{17};}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {  Surface{65};}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {  Surface{110};}

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {  Surface{155};}
