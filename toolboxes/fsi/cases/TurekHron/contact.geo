h=0.1;
h1=0.01;

Point(1) = {0., 0., 0., h};
Point(2) = {4, 0., 0., h};
Point(3) = {4, 2, 0., h};
Point(4) = {0., 2, 0., h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


Point(5) = {0.3, 1, 0., h1};
Point(6) = {0.3966,1 - 0.025,0.,h1};
Point(7) = {0.3966,1.025,0.,h1};
Point(8) = {0.2,1,0.,h1};

Circle(5) = {7, 5, 8};
Circle(6) = {8, 5, 6};

Point(9) = {2.5,1-0.025,0.,h1};
Point(10) = {2.5,1+0.025,0.,h1};
Point(11) = {2.5,1,0,h1};
Point(12) = {3,1.,0,h1*3.};

Line(7) = {6, 9};
Line(8) = {9, 11};
Line(9) = {11, 10};
Line(10) = {10, 7};
Line(11) = {11, 12};


// surface fluid
Line Loop(11) = {3, 4, 1, 2};
Line Loop(12) = {5, 6, 7, 8, 9, 10,11,-11};
Plane Surface(1) = {11, 12};

// surface structure
Circle(14) = {7, 5, 6};
Line Loop(15) = {10, 14, 7, 8, 9};
Plane Surface(2) = {15};

Physical Line("FEELPP_GMSH_PHYSICALNAME_IGNORED") = {11};

Physical Line("fluid-inlet") = {4};
Physical Line("fluid-wall") = {1,3};
Physical Line("fluid-cylinder") = {5,6};
Physical Line("fluid-outlet") = {2};
Physical Line("contact") = {7};
Physical Line("fsi-wall") = {8,9,10};
Physical Line("solid-fixed") = {14};

Physical Surface("Fluid") = {1};
Physical Surface("Solid") = {2};
