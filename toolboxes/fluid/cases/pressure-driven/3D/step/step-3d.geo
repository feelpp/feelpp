SetFactory("OpenCASCADE");

h=0.04;

hin=1;
hout=2;
w=2;
lstep=1;
lfloor=5;

Point(1) = {-lstep, hout-hin, 0, h};
Point(2) = {0, hout-hin, 0, h*0.1};
Point(3) = {0, 0, 0, h*0.5};
Point(4) = {lfloor, 0, 0, h};
Point(5) = {lfloor, hout, 0, h};
Point(6) = {-lstep, hout, 0, h};

Point(7) = {-lstep, hout-hin, w, h};
Point(8) = {0, hout-hin, w, h*0.1};
Point(9) = {0, 0, w, h*0.5};
Point(10) = {lfloor, 0, w, h};
Point(11) = {lfloor, hout, w, h};
Point(12) = {-lstep, hout, w, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 7};

Line(13) = {1, 7};
Line(14) = {2, 8};
Line(15) = {3, 9};
Line(16) = {4, 10};
Line(17) = {5, 11};
Line(18) = {6, 12};

// lateral walls
Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
Line Loop(2) = {7, 8, 9, 10, 11, 12};
Plane Surface(2) = {2};

// top and bottom walls
Line Loop(3) = {1, 14, -7, -13};
Plane Surface(3) = {3};
Line Loop(4) = {2, 15, -8, -14};
Plane Surface(4) = {4};
Line Loop(5) = {3, 16, -9, -15};
Plane Surface(5) = {5};
Line Loop(6) = {5, 18, -11, -17};
Plane Surface(6) = {6};

// inlet
Line Loop(7) = {-6, 18, 12, -13};
Plane Surface(7) = {7};

// outlet
Line Loop(9) = {10, -17, -4, 16};
Plane Surface(8) = {9};

Surface Loop(1) = {7, 1, 3, 4, 5, 8, 2, 6};
Volume(1) = {1};

// physical labels
Physical Surface("wall") = {1, 2, 3, 4, 5, 6};
Physical Surface("inlet") = {7};
Physical Surface("outlet") = {8};
Physical Volume("Omega") = {1};