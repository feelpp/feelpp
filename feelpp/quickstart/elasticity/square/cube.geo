SetFactory("OpenCASCADE");
Nc = 40;
Box(1) = {0, 0, 0, 1, 1, 1};

Transfinite Curve{1} = Nc;
Transfinite Curve{2} = Nc;
Transfinite Curve{3} = Nc;
Transfinite Curve{4} = Nc;
Transfinite Curve{5} = Nc;
Transfinite Curve{6} = Nc;
Transfinite Curve{7} = Nc;
Transfinite Curve{8} = Nc;
Transfinite Curve{9} = Nc;
Transfinite Curve{10} = Nc;
Transfinite Curve{11} = Nc;
Transfinite Curve{12} = Nc;

Physical Volume("body", 13) = {1};
Physical Surface("contact", 14) = {1};
Physical Surface("dirichlet", 15) = {2};
Physical Surface("wall", 16) = {4, 6, 5, 3};