
h=0.01;
hInf=h*10;//h*5;

r1=0.1;
r2=0.15;
height=0.10;
r_inf=2;//1;
r_ext=1;//0.7;

Point(1) = {0,0,0, h};
Point(2) = {r1,0,-height/2., h};
Point(3) = {r2,0,-height/2., h};
Point(4) = {r2,0, height/2., h};
Point(5) = {r1,0, height/2., h};

Point(6) = {r_ext,0,0, hInf};
Point(7) = {0,r_ext,0, hInf};
Point(8) = {0,-r_ext,0, hInf};
Point(9) = {-r_ext,0,0, hInf};
Point(10) = {0,0,r_ext, hInf};
Point(11) = {0,0,-r_ext, hInf};

// Define torus section and volume extruded
Line(1) = {5, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 5};
Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{1}; 
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2.} {
  Surface{26}; 
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2.} {
  Surface{48}; 
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2.} {
  Surface{70}; 
}

// Define Air 
Circle(87) = {6, 1, 7};
Circle(88) = {6, 1, 11};
Circle(89) = {7, 1, 11};
Circle(90) = {11, 1, 8};
Circle(91) = {8, 1, 10};
Circle(92) = {10, 1, 7};
Circle(93) = {6, 1, 10};
Circle(94) = {6, 1, 8};
Circle(95) = {8, 1, 9};
Circle(96) = {9, 1, 7};
Circle(97) = {11, 1, 9};
Circle(98) = {9, 1, 10};
Curve Loop(2) = {90, 95, -97};
Surface(92) = {2};
Curve Loop(3) = {97, 96, 89};
Surface(93) = {3};
Curve Loop(4) = {87, 89, -88};
Surface(94) = {4};
Curve Loop(5) = {90, -94, 88};
Surface(95) = {5};
Curve Loop(6) = {95, 98, -91};
Surface(96) = {6};
Curve Loop(7) = {98, 92, -96};
Surface(97) = {7};
Curve Loop(8) = {87, -92, -93};
Surface(98) = {8};
Curve Loop(9) = {94, 91, -93};
Surface(99) = {9};
Surface Loop(1) = {93, 92, 95, 99, 96, 97, 98, 94};
Surface Loop(2) = {35, 39, 43, 47, 57, 61, 65, 69, 21, 25, 13, 17, 87, 91, 79, 83};
Volume(5) = {1, 2};

// Markers
Physical Volume("Coil") = {1,2,3,4}; // Torus
Physical Volume("Omega") = {5}; // Infinite
Physical Surface("BoundaryInfinite") = {94, 98, 97, 93, 96, 99, 92, 95};
