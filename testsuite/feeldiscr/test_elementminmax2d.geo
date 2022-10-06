SetFactory("OpenCASCADE");
h = 0.1;
Rectangle(1) = {-0.5, -0.5, 0, 1, 1, 0};
Disk(2) = {0, 0, 0, 0.25, 0.25};

S[]=BooleanFragments{ Surface{1}; Delete; }{Surface{2}; Delete;};
Characteristic Length{ PointsOf{ Surface{ : }; } } = h;

Physical Curve("OuterBoundary") = {4, 3, 2, 1};
Physical Curve("InnerBoundary") = {5};
Physical Surface("Inner") = {2};
Physical Surface("Outer") = {3};
