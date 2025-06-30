h = 0.0002;
hBbox = h*10;
radiusSphere = 0.0001;
radiusConeLarge = 0.001;
lengthCone = 0.0075;
lengthHead = 0.0005;

// demi sphere
Point(1) = {-0.00775, 0, 0, h};
Point(2) = {-0.00775-radiusSphere, 0, 0, h};
Point(3) = {-0.00775, -radiusSphere, 0, h};
Point(4) = {-0.00775, radiusSphere, 0, h};
Point(5) = {-0.00775, 0, -radiusSphere , h};
Point(6) = {-0.00775, 0, radiusSphere , h};

// cone
Point(7) = {-0.00775+lengthCone, 0, 0, h};
Point(8) = {-0.00775+lengthCone, -radiusConeLarge, 0, h};
Point(9) = {-0.00775+lengthCone, radiusConeLarge, 0, h};
Point(10) = {-0.00775+lengthCone, 0, -radiusConeLarge , h};
Point(11) = {-0.00775+lengthCone, 0, radiusConeLarge , h};

// cylinder
Point(12) = {-0.00775+lengthCone+lengthHead, 0, 0, h};
Point(13) = {-0.00775+lengthCone+lengthHead, -radiusConeLarge, 0, h};
Point(14) = {-0.00775+lengthCone+lengthHead, radiusConeLarge, 0, h};
Point(15) = {-0.00775+lengthCone+lengthHead, 0, -radiusConeLarge , h};
Point(16) = {-0.00775+lengthCone+lengthHead, 0, radiusConeLarge , h};

// bbox
//Box(4) = {-0.02, -0.01, -0.01, 0.04, 0.02, 0.02};//+//+
bboxMinX = -0.02;
bboxMinY = -0.01;
bboxMinZ = -0.01;
bboxMaxX = 0.04;
bboxMaxY = 0.01;
bboxMaxZ = 0.01;
Point(17) = { bboxMinX, bboxMinY, bboxMinZ, hBbox};
Point(18) = { bboxMinX, bboxMinY, bboxMaxZ, hBbox};
Point(19) = { bboxMinX, bboxMaxY, bboxMinZ, hBbox};
Point(20) = { bboxMinX, bboxMaxY, bboxMaxZ, hBbox};
Point(21) = { bboxMaxX, bboxMinY, bboxMinZ, hBbox};
Point(22) = { bboxMaxX, bboxMinY, bboxMaxZ, hBbox};
Point(23) = { bboxMaxX, bboxMaxY, bboxMinZ, hBbox};
Point(24) = { bboxMaxX, bboxMaxY, bboxMaxZ, hBbox};


Circle(4) = {6, 1, 2};
Circle(5) = {6, 1, 3};
Circle(6) = {5, 1, 2};
Circle(7) = {5, 1, 3};
Circle(8) = {5, 1, 4};
Circle(9) = {6, 1, 4};
Circle(10) = {4, 1, 2};
Circle(11) = {2, 1, 3};

Curve Loop(1) = {9, 10, -4};
Surface(1) = {1};
Curve Loop(2) = {4, 11, -5};
Surface(2) = {2};
Curve Loop(3) = {11, -7, 6};
Surface(3) = {3};
Curve Loop(4) = {8, 10, -6};
Surface(4) = {4};



Circle(12) = {8, 7, 10};
Circle(13) = {10, 7, 9};
Circle(14) = {9, 7, 11};
Circle(15) = {11, 7, 8};
Circle(16) = {13, 12, 15};
Circle(17) = {15, 12, 14};
Circle(18) = {14, 12, 16};
Circle(19) = {16, 12, 13};
Line(20) = {6, 11};
Line(21) = {3, 8};
Line(22) = {5, 10};
Line(23) = {4, 9};
Line(24) = {11, 16};
Line(25) = {9, 14};
Line(26) = {8, 13};
Line(27) = {10, 15};
Curve Loop(5) = {14, 15, 12, 13};
Plane Surface(5) = {5};
Curve Loop(6) = {18, 19, 16, 17};
Plane Surface(6) = {6};
Curve Loop(7) = {20, 15, -21, -5};
Surface(7) = {7};
Curve Loop(8) = {21, 12, -22, 7};
Surface(8) = {8};
Curve Loop(9) = {8, 23, -13, -22};
Surface(9) = {9};
Curve Loop(10) = {23, 14, -20, 9};
Surface(10) = {10};
Curve Loop(11) = {13, 25, -17, -27};
Surface(11) = {11};
Curve Loop(12) = {25, 18, -24, -14};
Surface(12) = {12};
Curve Loop(13) = {27, -16, -26, 12};
Surface(13) = {13};
Curve Loop(14) = {19, -26, -15, 24};
Surface(14) = {14};
Surface Loop(1) = {9, 4, 1, 10, 7, 8, 3, 2, 5};
Volume(1) = {1};
Surface Loop(2) = {11, 12, 6, 14, 13, 5};
Volume(2) = {2};


Line(28) = {17, 19};
Line(29) = {19, 20};
Line(30) = {20, 18};
Line(31) = {18, 17};
Line(32) = {23, 24};
Line(33) = {24, 22};
Line(34) = {22, 21};
Line(35) = {21, 23};
Line(36) = {23, 19};
Line(38) = {21, 17};
Line(39) = {18, 22};
Line(40) = {20, 24};
Curve Loop(15) = {39, 34, 38, -31};
Plane Surface(15) = {15};
Curve Loop(16) = {38, 28, -36, -35};
Plane Surface(16) = {16};
Curve Loop(17) = {29, 30, 31, 28};
Plane Surface(17) = {17};
Curve Loop(18) = {34, 35, 32, 33};
Plane Surface(18) = {18};
Curve Loop(19) = {39, -33, -40, 30};
Plane Surface(19) = {19};
Curve Loop(20) = {36, 29, 40, -32};
Plane Surface(20) = {20};
Surface Loop(3) = {16, 15, 19, 18, 20, 17};
Surface Loop(4) = {9, 4, 1, 10, 7, 8, 3, 2, 6, 12, 11, 13, 14};
Volume(3) = {3, 4};

/*
Characteristic Length{ PointsOf{ Volume{5}; } } = h;
Characteristic Length{ PointsOf{ Volume{13}; } } = h;
*/

Physical Volume("Fluid") = {3};
Physical Volume("Solid") = {1};
Physical Volume("Head") = {2};
Physical Surface("fsi-wall") = {6, 12, 10, 7, 8, 9, 13, 14, 11, 2, 3, 4, 1};
Physical Surface("magneto") = {5};
Physical Surface("fluid-inlet") = {17};
Physical Surface("fluid-outlet") = {18};
Physical Surface("fluid-wall") = {20, 16, 15, 19};
