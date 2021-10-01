SetFactory("OpenCASCADE");

h = 1; // or whatever

Box(1) = {-5,0,-5,10,10,10};
Box(2) = {-0.05,8,-0.05,0.1,0.2,0.1};


V[]=BooleanFragments{Volume{1}; Delete; }{ Volume{2}; Delete;} ;

// Assign a mesh size to all the points of all the volumes:
MeshSize{ PointsOf{ Volume{3}; } } = h;

// Override this constraint on the points of the five spheres:
MeshSize{ PointsOf{ Volume{2}; } } = h/100;

Physical Surface("Wall_Body")={7:12};
Physical Surface("Wall")={13:18};
Physical Volume("Body")={2};
Physical Volume("Fluid")={3};