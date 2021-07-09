h=5;

Point(1) = {-1,-1,-1,h};
Point(2) = {1,-1,-1,h};
Point(3) = {-1,1,-1,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};

Line Loop(1) = {1,2,3};

Plane Surface(1) = {1};

Physical Line("Dirichlet") = {1,2,3};
Physical Surface("omega") = {1};

