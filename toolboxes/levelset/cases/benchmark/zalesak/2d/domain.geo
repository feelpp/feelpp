// In this file we describe the full domain, which is a square

// mesh size
h=0.05;

// square bounds
xmin = 0.0;
xmax = 1.0;
ymin = 0.0;
ymax = 1.0;

// square vertices
Point(1) = {xmin,ymin,0.0,h};
Point(2) = {xmax,ymin,0.0,h};
Point(3) = {xmax,ymax,0.0,h};
Point(4) = {xmin,ymax,0.0,h};

// square edges
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// square perimeter
Line Loop(5) = {1,2,3,4};

// square surface
Plane Surface(1) = {5};

// physical edges
Physical Line("Bottom") = {1};
Physical Line("Right") = {2};
Physical Line("Top") = {3};
Physical Line("Left") = {4};


// physical surface
Physical Surface("Omega") = {1};