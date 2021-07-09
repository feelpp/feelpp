h = 0.1;
L=1;
H=0.0635;
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
Point(1) = {0,0,0.0,h};
Point(2) = {L,0,0.0,h};
Point(3) = {L,H,0.0,h};
Point(4) = {0,H,0.0,h};
Line(1) = {4,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};


//Transfinite Line {1,3} = 41 Using Bump 0.1;

//Transfinite Line {1,3} = 20 Using Bump 0.1;
//Transfinite Line {2, 4} = 60 Using Progression 1;

/*
Transfinite Line {1,3} = 80 Using Bump 0.1;
Transfinite Line {2, 4} = 240 Using Progression 1;
Transfinite Surface {6};// Alternated;
*/
Physical Line("Gamma4") = {1};
Physical Line("Gamma3") = {2};
Physical Line("Gamma2") = {3};
Physical Line("Gamma1") = {4};
Physical Surface("Omega") = {6};
Physical Point("Corners") = {1,4};//{1,2,3,4};
