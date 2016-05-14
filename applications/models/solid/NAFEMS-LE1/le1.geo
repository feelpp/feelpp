h=0.1;

Point(1) = {0,0,0,h};
Point(2) = {0,1,0,h};
Point(3) = {2,0,0,h};
Point(4) = {3.25,0,0,h};
Point(5) = {0,2.75,0,h};

Ellipsis(1) = {2,1,2,3};
Line(2) = {3,4};
Ellipsis(3) = {4,1,4,5};
Line(4) = {5,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(1) = {5};

Physical Line( "AB" ) = {4};
Physical Line( "BC" ) = {3};
Physical Line( "CD" ) = {2};
Physical Line( "DA" ) = {1};

Physical Surface( "Omega" ) = {1};
