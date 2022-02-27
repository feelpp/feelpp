
h=0.01;
hFine=h/10;
/*
scale=0.0127;
//L=2;
L1=110*scale;
L2=50*scale;
InletHeight=8*scale;
StepHeight=1*scale;
*/
L1=0.3048;
L2=1.3335;
StepHeight=0.0381;
InletHeight=0.0762;

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

Point(10) = {L2, StepHeight+InletHeight/2, 0, hFine};
Point(11) = {-L1, StepHeight+InletHeight/2, 0, h};
Point(12) = {0, StepHeight+InletHeight/2, 0, hFine};

Line(1) = {8, 9};
Line(2) = {9, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 10};
Line(6) = {10, 5};
Line(7) = {6, 5};
Line(8) = {6, 7};
Line(9) = {7, 11};
Line(10) = {11, 8};
Line(11) = {6, 12};
Line(12) = {12, 9};
Line(13) = {9, 4};
Line(14) = {11, 12};
Line(15) = {12, 10};
Line Loop(16) = {9, 14, -11, 8};
Plane Surface(17) = {16};
Line Loop(18) = {11, 15, 6, -7};
Plane Surface(19) = {18};
Line Loop(20) = {15, -5, -13, -12};
Plane Surface(21) = {20};
Line Loop(22) = {12, -1, -10, 14};
Plane Surface(23) = {22};
Line Loop(24) = {13, -4, -3, -2};
Plane Surface(25) = {24};

Physical Line("inlet") = {9,10};
Physical Line("wall_horizontal") = {1,3,7,8};
Physical Line("wall_vertical") = {2};
Physical Line("outlet") = {4,5,6};
Physical Surface("fluid") = {17,23,19,21,25};

Transfinite Line {1,14,8} = 240 Using Progression 1;
Transfinite Line {3,13,15,7} = 150 Using Progression 1.02;
Transfinite Line {9,11,6,10,12,5} = 40 Using Bump 0.05;// 0.1;
Transfinite Line {2,4} = 40 Using Bump 0.05;//0.1;
Transfinite Surface {17} Right;// Alternated;
Transfinite Surface {23} Right;
Transfinite Surface {19} Right;
Transfinite Surface {21} Right;
Transfinite Surface {25} Right;