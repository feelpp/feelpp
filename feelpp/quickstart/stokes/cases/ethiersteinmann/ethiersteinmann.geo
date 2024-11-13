//+
SetFactory("OpenCASCADE");
//+
h=0.2;
Box(1) = {0, 0, 0, 1, 1, 1};
Characteristic Length{ PointsOf{ Volume{:}; } } = h;
//+
Physical Surface("Dirichlet") = {6, 5};
//+
Physical Surface("Neumann") = {1, 2, 4, 3};
//+
Physical Volume("Omega") = {1};