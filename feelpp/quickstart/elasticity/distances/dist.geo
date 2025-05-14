SetFactory("OpenCASCADE");
h = 0.1;
Box(1) = {0, 0, 0, 1, 1, 1};

Characteristic Length { PointsOf{ Volume{1}; } } = h;

Physical Surface("boundary", 13) = {6, 5, 4, 1, 2, 3};
Physical Volume("volume", 14) = {1};

