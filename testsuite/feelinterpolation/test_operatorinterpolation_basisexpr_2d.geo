h = 0.05;//1e-3;
R=0.5;//0.001;
ratio=1;
S=2*R*Sqrt(1/ratio);

Point(0) = {0,0,0,h};


Point(1) = {Sqrt(ratio)*R, 0, 0, h};
Point(2) = {0, R/Sqrt(ratio), 0, h};
Point(3) = {-Sqrt(ratio)*R, 0, 0, h};
Point(4) = {0, -R/Sqrt(ratio), 0, h};

Point(5) = {S, 0, 0, 0.5*h};
Point(6) = {-S, 0, 0, 0.5*h};

Circle(1) = {1, 0, 2};
Circle(2) = {2, 0, 3};
Circle(3) = {3, 0, 4};
Circle(4) = {4, 0, 1};

Circle(5) = {5, 0, 6};
Circle(6) = {6, 0, 5};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {5,6};
Plane Surface(2) = {2,1};

Physical Curve("Interface") = {1,2,3,4};
Physical Surface("Omega2") = {1};
Physical Curve("ExteriorBoundary") = {5,6};
Physical Surface("Omega1") = {2};


