lc = 1;
Point(1) = {-1,-1,-1,lc};
Point(2) = { 1,-1,-1,lc};
Point(3) = { 1, 1,-1,lc};
Point(4) = {-1, 1,-1,lc};
Line(1) = {4,3};
Line(2) = {3,2};
Line(3) = {2,1};
Line(4) = {1,4};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
tmp[] += Extrude {0,0.0,2} {
    Surface{6};
};
Physical Surface("Border") = {27,19,23,15,28,6};
Physical Volume("Omega") = tmp[1];
