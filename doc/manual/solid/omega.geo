h = 0.1;
L=10;
H=1;
Point(1) = {0, -H/2, 0, h};
Point(2) = {0, H/2, 0, h};
Point(3) = {L, H/2, 0, h};
Point(4) = {L, -H/2, 0, h};
Point(5) = {0,0.,0,h};
Point(6) = {0.1,0.1,0.,h};

Line(1) = {2, 5};

Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line(11) = {5,1};

Curve Loop(1) = {1, 11, 2, 3, 4};
Plane Surface(1) = {1};
// Include points in the surface
//Point{6} In Surface{1};

//Physical Curve("Wall", 8) = {1, 11,2, 3, 4};
Physical Surface("Caoutchouc", 9) = {1};
//Physical Curve("Upper", 10) = {4};
//Physical Curve("Lower", 11) = {2};
//Physical Point("a0") = {5};
Physical Point("a1") = {6};
Physical Point("Pinned") = {5}; // Specify the pinned end
Physical Curve("FreeEnd") = {3}; // Specify the free end for force