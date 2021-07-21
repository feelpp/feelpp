
h=0.01;
hFine=h/10;
scale=0.0127;
//L=2;
L1=110*scale;
L2=50*scale;
InletHeight=8*scale;
StepHeight=1*scale;

OutletHeight=StepHeight+InletHeight;
//Point(1) = {-1, -1, 0, h};
Point(2) = {0, 0, 0, hFine};
Point(3) = {L2, 0, 0, hFine};
Point(4) = {L2, StepHeight, 0, hFine};
Point(5) = {L2, OutletHeight, 0, h};
Point(6) = {0, OutletHeight, 0, h};
Point(7) = {-L1, OutletHeight, 0, h};
Point(8) = {-L1, StepHeight, 0, h};
Point(9) = {0, StepHeight, 0, hFine};
Line(1) = {8, 9};
Line(2) = {9, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {9, 4};

Line Loop(10) = {8, 1, 9, 5, 6, 7};
Plane Surface(11) = {10};
Line Loop(12) = {9, -4, -3, -2};
Plane Surface(13) = {12};

Physical Line("inlet") = {8};
Physical Line("wall") = {6, 7, 1, 2, 3};
Physical Line("outlet") = {5, 4};
Physical Surface("fluid") = {11,13};



/*
// boundary layer
Field[1] = BoundaryLayer;
Field[1].EdgesList = {6,7};
Field[1].hfar = 0.01;
Field[1].hwall_n = 0.001;
Field[1].thickness = 0.02;
Field[1].ratio = 1.1;
Field[1].IntersectMetrics = 1;
Field[1].Quads = 0;//1;
BoundaryLayer Field = 1;
Transfinite Line{8} = 30 Using Bump 0.02;
*/
