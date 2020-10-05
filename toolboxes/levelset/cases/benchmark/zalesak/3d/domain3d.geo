// In this file we describe the full domain, which is a unit cube

// mesh size
h=0.05;

// cube bounds
xmin = 0.0;
xmax = 1.0;
ymin = 0.0;
ymax = 1.0;
zmin = 0.0;
zmax = 1.0;

// square vertices (bottom of the cube)
Point(1) = {xmin,ymin,zmin,h};
Point(2) = {xmax,ymin,zmin,h};
Point(3) = {xmax,ymax,zmin,h};
Point(4) = {xmin,ymax,zmin,h};
// square vertices (top of the cube)
Point(5) = {xmin,ymin,zmax,h};
Point(6) = {xmax,ymin,zmax,h};
Point(7) = {xmax,ymax,zmax,h};
Point(8) = {xmin,ymax,zmax,h};

// cube edges (bottom square)
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
// cube edges (vertical edges)
Line(5) = {1,5};
Line(6) = {2,6};
Line(7) = {3,7};
Line(8) = {4,8};
// cube edges (top square)
Line(9) = {5,6};
Line(10) = {6,7};
Line(11) = {7,8};
Line(12) = {8,5};

// squares (cube faces) perimeters
// bottom
Line Loop(13) = {1,2,3,4};
// front
Line Loop(14) = {1,6,-9,-5};
// right
Line Loop(15) = {2,7,-10,-6};
// back
Line Loop(16) = {3,8,-11,-7};
// left
Line Loop(17) = {4,5,-12,-8};
// top
Line Loop(18) = {9,10,11,12};

// squares (cube faces) surfaces
Plane Surface(1) = {13};
Plane Surface(2) = {14};
Plane Surface(3) = {15};
Plane Surface(4) = {16};
Plane Surface(5) = {17};
Plane Surface(6) = {18};

// cube surface
Surface Loop(7) = {1,2,3,4,5,6};

// cube volume
Volume(1) = {7};

// physical edges
//Physical Line("Test") = {1};

// physical surfaces
Physical Surface("Bottom") = {1};
Physical Surface("Front") = {2};
Physical Surface("Right") = {3};
Physical Surface("Back") = {4};
Physical Surface("Left") = {5};
Physical Surface("Top") = {6};

// Physical volume
Physical Volume("Omega") = {1};