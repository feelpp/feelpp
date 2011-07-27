h=0.25;
Point(1) = {0, 0, 0, h};
Point(2) = {0, 1, 0, h};
Point(3) = {0, 2, 0, h};
Point(4) = {2, 0, 0, h};
Point(5) = {2, 1, 0, h};
Point(6) = {2, 2, 0, h};
Point(7) = {3, 0, 0, h};
Point(8) = {4, 0, 0, h};
Point(9) = {4, 2, 0, h};
Point(10) = {2.707106781, 0.707106781, 0, h};

Line(1) = {2, 3};
Line(2) = {3, 6};
Line(3) = {9, 6};
Line(4) = {9, 8};
Line(5) = {8, 7};
Line(6) = {10, 9};
Line(7) = {6, 5};
Line(8) = {5, 2};
Circle(9) = {5, 4, 10};
Circle(10) = {7, 4, 10};
//Physical Line(2) = {1}
//Physical Line(1) = {2,3,4,5,8,9,10}

Physical Line(2) = {1,2,3,4,5,8,9,10}

//Physical Line(8) = {1,2,3,4,5,9,10,8};
//Physical Line(8) = {1};
Line Loop(1) = {1, 8, 7, 2};
Plane Surface(1) = {1};
Line Loop(2) = {7, 9, 6, 3};
Plane Surface(2) = {2};
Line Loop(3) = {4, 6, 10, 5};
Plane Surface(3) = {3};
Physical Surface(9) = {-1};
Physical Surface(10) = {2};
Physical Surface(11) = {-3};


