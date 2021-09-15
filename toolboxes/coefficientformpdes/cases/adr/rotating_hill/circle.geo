h=0.05;

SetFactory("OpenCASCADE");
Characteristic Length {1} = 0.01;

Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
Curve Loop(1) = {1};
Plane Surface(1) = {1};
Characteristic Length{ PointsOf{ Surface{1}; } } = h;

Physical Line("Gamma") = {1};
Physical Surface("Omega") = {1};
