Nc = 5;

Point(1) = {0, 0, 0, 1};
Point(2) = {0, 1, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {1, 0, 0, 1};

Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};

Curve Loop(1) = {4, 1, 2, 3};
Surface(1) = {1};

Transfinite Curve{1} = Nc;
Transfinite Curve{2} = Nc;
Transfinite Curve{3} = Nc;
Transfinite Curve{4} = Nc;

Physical Curve("contact", 5) = {4};
Physical Curve("dirichlet", 6) = {2};
Physical Curve("wall", 7) = {3, 1};
Physical Surface("body", 8) = {1};