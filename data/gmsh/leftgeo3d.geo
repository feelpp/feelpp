h = 0.1;

Point(1) = {0,0,0,h};
Point(2) = {1,0,0,h};
Point(3) = {0,1,0,h};
Point(4) = {0,0,1,h};
Point(5) = {-1,0,0,h};
Point(6) = {0,-1,0,h};
Point(7) = {0,0,-1,h};

Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,6};
Circle(4) = {6,1,2};
Circle(5) = {2,1,7};
Circle(6) = {7,1,5};
Circle(7) = {5,1,4};
Circle(8) = {4,1,2};
Circle(9) = {6,1,7};
Circle(10) = {7,1,3};
Circle(11) = {3,1,4};
Circle(12) ={4,1,6};

Physical Line(10) = {1};

Line Loop(1) = {1,11,8};
Line Loop(2) = {2,7,-11};
Line Loop(3) = {3,-12,-7};
Line Loop(4) = {4,-8,12};
Line Loop(5) = {5,10,-1};
Line Loop(6) = {-2,-10,6};
Line Loop(7) = {-3,-6,-9};
Line Loop(8) = {-4,9,-5};

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};
Ruled Surface(7) = {7};
Ruled Surface(8) = {8};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};
Physical Surface(8) = {8};

Surface Loop (1) = {1,2,3,4,5,6,7,8};

Volume (1) = {1};

Physical Volume("Mat1") = {1};
