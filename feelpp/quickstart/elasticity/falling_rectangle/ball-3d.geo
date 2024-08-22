SetFactory("OpenCASCADE");

h = 1;

Sphere(1) = {0, 0, 24, 20};

Point(100) = {0, 0, 24, h};

Point{100} In Volume{1};

Physical Point("CM")={100};
Physical Surface("Wall") = {1};
Physical Volume("Caoutchouc") = {1};

hwall = 1;

Box(2) = {-25, -25, -1, 50, 50, 1};

Physical Surface("Obs1") = {2,3,4,5,6,7};
Physical Volume("obstacle") = {2};

Characteristic Length { PointsOf{ Volume{1}; } } = h;
Characteristic Length { PointsOf{ Volume{2}; } } = hwall;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;