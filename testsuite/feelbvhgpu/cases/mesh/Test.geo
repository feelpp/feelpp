//+
SetFactory("OpenCASCADE");
//+
Box(1) = {-2, -2, -2, 4, 4, 4};
//+
Box(2) = {-1, -1, -1, 2, 2, 2};
//+
Box(3) = {-3.5, -0.5, -0.5, 0.5, 0.5, 0.5};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Physical Volume("Material", 25) = {1};
//+
Physical Surface("CavityBottom", 26) = {8};
//+
Physical Surface("CavityTop", 27) = {10};
//+
Physical Surface("CavitySides", 28) = {7, 9, 12, 11};
//+
Physical Volume("Material2", 49) = {3};


//+
Physical Surface("Up3", 50) = {16};
//+
Physical Surface("Down3", 51) = {15};
//+
Physical Surface("Front3", 52) = {18};
//+
Physical Surface("Back3", 53) = {17};
//+
Physical Surface("Left3", 54) = {13};
//+
Physical Surface("Rigth3", 55) = {14};
