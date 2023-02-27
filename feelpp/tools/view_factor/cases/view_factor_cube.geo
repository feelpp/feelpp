// Gmsh project created on Wed Mar 23 14:24:58 2022

h=0.05;

SetFactory("OpenCASCADE");
//+
//Rectangle(1) = {0, 0, 0, 1, 0.5, 0};
//+
//Rectangle(2) = {0, 0, 1, 1, 0.5, 0};
//+
//Physical Surface("TopWall") = {2};
//+
//Physical Surface("BottomWall") = {1};
//+
//Box(1) = {0, 0, 0, 1, 1, 1};
Box(1) = {-1,-1,-0.8, 2, 2, 2};
//+
Physical Surface("TopWall") = {6};
//+
Physical Surface("Walls") = {1, 4, 2, 3};
//+
Physical Surface("BottomWall") = {5};
//+
Physical Volume("Air") = {1};

Characteristic Length {4, 8, 7, 3, 2, 6, 1, 5} = h;
