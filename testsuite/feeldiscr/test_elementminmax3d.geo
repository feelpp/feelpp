SetFactory("OpenCASCADE");
h=0.1;
Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};
Sphere(2) = {0, 0, 0, 0.25, -Pi/2, Pi/2, 2*Pi};

S[]=BooleanFragments{ Volume{1}; Delete; }{Volume{2}; Delete;};
Characteristic Length{ PointsOf{ Volume{ : }; } } = h;

Physical Surface("OuterBoundary") = {3, 1, 4, 2,5,6};
Physical Surface("InnerBoundary") = {7};
Physical Volume("Outer") = {4};
Physical Volume("Inner") = {2};