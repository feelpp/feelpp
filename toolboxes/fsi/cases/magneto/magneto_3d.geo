SetFactory("OpenCASCADE");
h = 0.0005;

Sphere(1) = {-0.00775, 0, 0, 0.0001, -Pi/2, Pi/2, Pi};
Rotate {{0, 0, 1}, {-0.00775, 0, 0}, Pi/2} {
  Volume{1};
}

Cone(2) = {-0.00775, 0, 0, 0.0075, 0, 0, 0.0001, 0.0015};

BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; }

Cylinder(3) = {-0.00025, 0, 0, 0.0005, 0, 0, 0.0015, 2*Pi};

Characteristic Length{ PointsOf{ Volume{1}; } } = h;
Characteristic Length{ PointsOf{ Volume{3}; } } = h;

Box(4) = {-0.02, -0.01, -0.01, 0.04, 0.02, 0.02};//+//+
BooleanDifference{ Volume{4}; Delete; }{ Volume{1}; Volume{3}; Delete; }

Sphere(11) = {-0.00775, 0, 0, 0.0001, -Pi/2, Pi/2, Pi};
Rotate {{0, 0, 1}, {-0.00775, 0, 0}, Pi/2} {
  Volume{11};
}

Cone(12) = {-0.00775, 0, 0, 0.0075, 0, 0, 0.0001, 0.0015};

BooleanUnion{ Volume{11}; Delete; }{ Volume{12}; Delete; }

Cylinder(13) = {-0.00025, 0, 0, 0.0005, 0, 0, 0.0015, 2*Pi};

Characteristic Length{ PointsOf{ Volume{5}; } } = h;
Characteristic Length{ PointsOf{ Volume{13}; } } = h;

Physical Volume("Fluid", 32) = {4};
Physical Volume("Solid", 33) = {5};
Physical Volume("Head", 34) = {13};
Physical Surface("fsi-wall", 35) = {13, 12, 5, 1};
Physical Surface("fluid-wall", 36) = {8, 9, 10, 7};
Physical Surface("fluid-outlet", 37) = {11};
Physical Surface("fluid-inlet", 38) = {6};
Physical Surface("magneto", 39) = {16};
