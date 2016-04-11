/*Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthExtendFromBoundary=1;
Mesh.CharacteristicLengthFromPoints=1;
//Mesh.ElementOrder=3;
Mesh.SecondOrderIncomplete = 0;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1; // 3D mesh algorithm (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D)
Mesh.RecombinationAlgorithm=0;
*/
Mesh.OptimizeNetgen=1;

//lc=0.2;//M0
//lc=0.15;//M1
//lc=0.1;//M2
//lc=0.075;//M3
//lc=0.05;//M4
h=0.05;
lc=h;
Mesh.CharacteristicLengthMax=lc;
Mesh.CharacteristicLengthMin=lc;


R=0;//1.5;
r=0.5;
//Point(1) = {0,0,0,lc};
Point(2) = {R,0,0,lc};
Point(3) = {R-r,0,0,lc};
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{3}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{4}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{5}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{6}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{7}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{8}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{9}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{10}; }

long=5;
Extrude {0, long, 0} { Line{2,1,8,7,6,5,3,4}; }


ep=0.1;
Point(1001) = {R-r-ep,0,0,lc}; 
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{1001}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{1002}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{1003}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{1004}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{1005}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{1006}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{1007}; }
Extrude {{0,1,0}, {R,0,0}, Pi/4} { Point{1008}; }

Extrude {0, long, 0} { Line{42,41,48,47,46,45,43,44}; }


//inlet fluid
Line Loop(81) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(82) = {81};
//inlet wall
Line Loop(83) = {41, 42, 43, 44, 45, 46, 47, 48};
Plane Surface(84) = {81, 83};

//outlet fluid
Line Loop(85) = {33, 37, 29, 25, 21, 17, 13, 9};
Plane Surface(86) = {85};
//outlet wall
Line Loop(87) = {73, 77, 69, 65, 61, 57, 53, 49};
Plane Surface(88) = {85,87};

// fluid volume
Surface Loop(89) = {86, 36, 82, 16, 12, 20, 24, 28, 32, 40};
Volume(90) = {89};
// wall volume
Surface Loop(91) = {72, 84, 60, 56, 52, 76, 80, 88, 64, 68, 32, 28, 24, 20, 16, 12, 36, 40};
Volume(92) = {91};

Physical Surface("inletBlood") = {82};
Physical Surface("outletBlood") = {86};
Physical Surface("fsiWall") = {12, 16, 20, 24, 28, 40, 32, 36};
Physical Surface("inletRing") = {84};
Physical Surface("outletRing") = {88};
Physical Surface("exterior") = {52, 76, 80, 72, 68, 64, 60, 56};

Physical Volume("ArterialWall")={92};
Physical Volume("Blood")={90};
