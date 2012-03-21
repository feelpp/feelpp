hs=0.24;
Mesh.ElementOrder = 4;
Point(1) = {0, 0, 0, hs};
Point(2) = {5, 0, 0, hs};
Point(3) = {0, 5, 0, hs};
Point(4) = {-5, 0, 0, hs};
Point(5) = {0, -5, 0, hs};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Surface("Bord C") = {1};

