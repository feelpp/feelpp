c=1;
hmax =0.0008;
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
Line(4) = {8, 7};
Line(5) = {7, 5};


Line(26) = {27, 28};
Line(27) = {29, 30};
Line(28) = {31, 32};
Line(29) = {33, 34};
Line(30) = {35, 36};
Line(31) = {37, 38};

Line Loop(32) = {4, 5, 1, -2};
Plane Surface(33) = {32};
