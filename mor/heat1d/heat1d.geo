// -*- mode: c++ -*-
h=0.1;
Point(1) = {-0.5,0,0,h};
Point(2) = {-0.1,0,0,h};
Point(3) = {0,0,0,h};
Point(4) = {0.1,0,0.0,h};
Point(5) = {0.5,0,0.0,h};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line Loop(5) = {1,2,3,4};

Physical Line("k1_1") = {1};
Physical Line("k1_2") = {2};
Physical Line("k2_1") = {3};
Physical Line("k2_2") = {4};

Physical Point("left") = {1};
Physical Point("right") = {5};
