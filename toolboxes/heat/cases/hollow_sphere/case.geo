
h = 0.05;

SetFactory("OpenCASCADE");

Sphere(1) = {0, 0, 0, 0.3, -Pi/12, Pi/12, Pi/6.};

Sphere(2) = {0, 0, 0, 0.392, -Pi/12, Pi/12, Pi/6};
Delete {
  Volume{1,2}; 
}

Recursive Delete {
  Surface{2,7,10,5,4,9,8,3}; 
}

Line(14) = {2, 8};
Line(15) = {7, 1};
Line(16) = {4, 10};
Line(17) = {3, 9};
Curve Loop(7) = {15, -1, 14, 10};
Surface(7) = {7};
Curve Loop(9) = {4, 14, -13, -16};
Surface(8) = {9};
Curve Loop(11) = {12, -17, -3, 16};
Surface(9) = {11};
Curve Loop(13) = {17, 11, 15, -2};
Surface(10) = {13};
Surface Loop(3) = {7, 6, 8, 9, 1, 10};
Volume(1) = {3};

Physical Volume("Omega") = {1};
Physical Surface("internal_surface") = {1};
Physical Surface("external_surface") = {6};
Physical Surface("others_surface") = {7, 8, 9, 10};

Characteristic Length {3, 1, 9, 7, 2, 8, 4, 10} = h;