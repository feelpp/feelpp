h=0.1;
hBody=h/3;
L=1;
Point(1) = {2.0000000000000001e-01,0.0000000000000000e+00,0.0000000000000000e+00,hBody};
Point(2) = {0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,hBody};
Point(3) = {-2.0000000000000001e-01,0.0000000000000000e+00,0.0000000000000000e+00,hBody};
Point(4) = {-1.0000000000000000e+00,-L,0.0000000000000000e+00,h};
Point(5) = {1.0000000000000000e+00,-L,0.0000000000000000e+00,h};
Point(6) = {1.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,h};
Point(7) = {-1.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,h};

Circle(1) = { 1,2,3};
Circle(2) = { 3,2,1};
Line(3) = { 4,5};
Line(4) = { 5,6};
Line(5) = { 6,7};
Line(6) = { 7,4};
Curve Loop(1) = {1,2};
Curve Loop(2) = {3,4,5,6};
Plane Surface(3) = {2,1};

Plane Surface(7) = {1};


Physical Curve("wall_bottom") = {3};
Physical Curve("wall_top") = {5};
Physical Curve("wall_left") = {6};
Physical Curve("wall_right") = {4};
Physical Curve("wall_body") = {1,2};
Physical Surface("Omega_fluid") = {3};
Physical Surface("Omega_body") = {7};
