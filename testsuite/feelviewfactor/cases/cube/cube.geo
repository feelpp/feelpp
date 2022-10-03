SetFactory("OpenCASCADE");
//+
hs = 0.1;

Box(1) = {-0.2, -0.2, -0.2, 1.4, 1.4, 1.4};
Box(2) = {0, 0, 0, 1, 1, 1};

S[]=BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete;};//+

Characteristic Length{ PointsOf{ Volume{ : }; } } = hs;
bdy[] = CombinedBoundary { Volume{1}; };

For ii In { 0 : (#bdy[]-1) }
        Printf("boundary number %g = %g", ii, Abs(bdy[ii]));
        Physical Surface(Sprintf("Gamma_%g",ii)) = Abs(bdy[ii]);   
EndFor
Physical Volume("Omega")={1};