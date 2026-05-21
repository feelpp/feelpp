SetFactory("OpenCASCADE");

h = 1;
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
zmin = 0;
epaisseur = 1;

Point(1) = {xmin, ymin, zmin, h};
Point(2) = {xmax, ymin, zmin, h};
Point(3) = {xmax, ymax, zmin, h};
Point(4) = {xmin, ymax, zmin, h};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve {1,2,3,4} = 1 Using Progression 1;
Recombine Surface {1};

Extrude {0, 0, epaisseur} {
    Surface{1}; Layers{1}; Recombine;
}

Transfinite Curve {1,2,3,4} = 1 Using Progression 1;



Physical Volume("Omega") = {1};

Physical Surface("ForceApply") = {3};
Physical Surface("Dirichlet") = {5};


// Physical Point("DirichletPoints") = {1,4,5,8};


// Physical Surface("XPlus") = {3};
// Physical Surface("XMoins") = {5};

// Physical Surface("YPlus") = {4};
// Physical Surface("YMoins") = {2};

// Physical Surface("ZPlus") = {6};
// Physical Surface("ZMoins") = {1};
