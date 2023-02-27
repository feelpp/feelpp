//+
SetFactory("OpenCASCADE");
//+
Box(1) = {-2, -2, -2, 4, 4, 4};
//+
Box(2) = {-1, -1, -1, 2, 2, 2};
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
