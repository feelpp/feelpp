h=1;

Point(1) = {-1,-1,-1,h};
Point(2) = {1,-1,-1,h};
Point(3) = {-1,1,-1,h};
Point(4) = {-1,-1,1,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};
Line(4) = {4,1};
Line(5) = {4,2};
Line(6) = {4,3};

Line Loop(1) = {1,2,3};
Line Loop(2) = {1,-5,4};
Line Loop(3) = {2,-6,5};
Line Loop(4) = {3,-4,6};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Surface Loop(1) = {1,2,3,4};
Volume(1) = {1};

//Physical Surface("Neumann") = {1};
Physical Surface("Dirichlet") = {1,2,3,4};
Physical Volume("omega") = {1};

