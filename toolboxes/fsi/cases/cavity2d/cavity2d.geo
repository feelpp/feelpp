h=0.05;
h1=0.003;//0.01;//0.003;
Point(1) = {0.0000000000000000e+00,0.0000000000000000e+00,0.0, h1};
Point(2) = {1.0000000000000000e+00,0.0000000000000000e+00,0.0, h1};
Point(3) = {1.0000000000000000e+00,8.7500000000000000e-01,0.0, h};
Point(4) = {1.0000000000000000e+00,1.0000000000000000e+00,0.0, h};
Point(5) = {0.0000000000000000e+00,1.0000000000000000e+00,0.0, h};
Point(6) = {0.0000000000000000e+00,8.7500000000000000e-01,0.0, h};

Point(7) = {0.0000000000000000e+00,-0.002,0.0, h1};
Point(8) = {1.0000000000000000e+00,-0.002,0.0, h1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line Loop(1) = {1,2,3,4,5,6};
Plane Surface(1) = {1};

Line(7) = {1,7};
Line(8) = {7,8};
Line(9) = {8,2};
Line Loop(2) = {7,8,9,-1};
Plane Surface(2) = {2};

Physical Line("fluid-inlet") = {5};
Physical Line("fluid-noslip-wall") = {2,6};
Physical Line("fluid-slip-wall") = {4};
Physical Line("fluid-outlet") = {3};
Physical Line("solid-fixed-wall") = {7,9};
Physical Line("solid-free-wall") = {8};

Physical Line("fsi-wall") = {1};

Physical Surface("Fluid") = {1};
Physical Surface("Solid") = {2};
