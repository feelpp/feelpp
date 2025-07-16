SetFactory("OpenCASCADE");

h = 0.0001;

height = 0.0015;
length_head = 0.0005;
length_tail = 0.0075;

height_box = 0.02;
lenght_box = 0.04;

Point(1) = {-length_head/2, -height/2, 0, h};
Point(2) = {-length_head/2, height/2, 0, h};
Point(3) = {length_head/2, height/2, 0, h};
Point(4) = {length_head/2, -height/2, 0, h};

Point(5) = {-length_head/2 - length_tail, -0.0001, 0, h};  // Point de départ de la queue
Point(6) = {-length_head/2 - length_tail, 0.0001, 0, h};  // Point de fin de la queue
Point(7) = {-length_head/2 - length_tail + 0.00005, 0, 0, h};  // Point pour début de l'arrondi


Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};

Line(28) = {1, 5};
Line(29) = {2, 6};

Circle(30) = {6, 7, 5};  

h1 = 0.001;

Point(11) = {-lenght_box/2, -height_box/2, 0, h1};
Point(12) = {lenght_box/2., -height_box/2, 0, h1};
Point(13) = {-lenght_box/2, height_box/2, 0, h1};
Point(14) = {lenght_box/2., height_box/2, 0, h1};

Line(5) = {13, 11};
Line(6) = {11, 12};
Line(7) = {12, 14};
Line(8) = {14, 13};

Curve Loop(1) = {29, 30, -28, -1};
Plane Surface(1) = {1};
Curve Loop(2) = {4, 1, 2, 3};
Plane Surface(2) = {2};
Curve Loop(3) = {8, 5, 6, 7};
Curve Loop(4) = {29, 30, -28, 2, 3, 4};
Plane Surface(3) = {3, 4};

Physical Surface("Solid", 31) = {1};
Physical Surface("Head", 32) = {2};
Physical Surface("Fluid", 33) = {3};
Physical Curve("fluid-inlet", 34) = {5};
Physical Curve("fluid-wall", 35) = {8, 6};
Physical Curve("fluid-outlet", 36) = {7};
Physical Curve("magneto", 37) = {1};
Physical Curve("fsi-wall", 38) = {29, 30, 28, 2, 3, 4};
