h = 0.1;
//Point(1) = {1, 0, 0, h};
Point(1) = {0, 0, 0, h};
Point(2) = {2, 0, 0, h};
Point(3) = {2, 1, 0, h};
Point(4) = {0, 1, 0, h};

//Point(5) = {1, 0, 1, h};
Point(5) = {0, 0, 1, h};
Point(6) = {2, 0, 1, h};
Point(7) = {2, 1, 1, h};
Point(8) = {0, 1, 1, h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line(9)  = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Physical Line(10) = {4};

// bottom
Line Loop(21) = {-1,-4,-3,-2};
Plane Surface(31) = {21} ;
Physical Surface(1) = {31};

// top
Line Loop(22) = {5,6,7,8};
Plane Surface(32) = {22} ;
Physical Surface(2) = {32};

// left
Line Loop(23) = {1,10,-5,-9};
Plane Surface(33) = {23} ;
Physical Surface(3) = {33};

// right
Line Loop(24) = {12,-7,-11,3};
Plane Surface(34) = {24} ;
Physical Surface(4) = {34};

// front
Line Loop(25) = {2,11,-6,-10};
Plane Surface(35) = {25} ;
Physical Surface(5) = {35};

// back
Line Loop(26) = {9,-8,-12,4};
Plane Surface(36) = {26} ;
Physical Surface(6) = {36};

Surface Loop(41) = {31,32,33,34,35,36};
Volume(51) = {41};

Physical Volume("Mat1") = {51};
