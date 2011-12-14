h=0.08;
Point(1) = {-0.8,-0.4,0,h};
Point(2) = {0.8,-0.4,0,h};
Point(3) = {0.8,0.4,0,h};
Point(4) = {-0.8,0.4,0,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Line Loop(1) = {2,3,4,1};

Plane Surface(3) = {1};


Physical Line("dirichlet") = {4};

Physical Line(12) = {2};
Physical Line(13) = {1,3};
Physical Surface(1) = {1,2,3,4};


