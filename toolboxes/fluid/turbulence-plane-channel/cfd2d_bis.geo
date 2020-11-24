h = 0.1;
L=1;
H=0.0635;
Point(1) = {0,0,0.0,h};
Point(2) = {L,0,0.0,h};
Point(3) = {L,H/2,0.0,h};
Point(4) = {L,H,0.0,h};
Point(5) = {0,H,0.0,h};
Point(6) = {0,H/2,0.0,h};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {3,6};

Line Loop(8) = {5, -7, 3, 4};
Plane Surface(9) = {8};
Line Loop(10) = {7, 6, 1, 2};
Plane Surface(11) = {10};

Transfinite Line {2,3,5,6} = 80 Using Bump 0.1;
Transfinite Line {1, 4,7} = 240 Using Progression 1;
Transfinite Surface {9} Right;// Alternated;
Transfinite Surface {11} Right;// Alternated;


//Line Loop(5) = {1,2,3,4};
//Plane Surface(6) = {5};

/*
//Transfinite Line {1,3} = 41 Using Bump 0.1;

//Transfinite Line {1,3} = 20 Using Bump 0.1;
//Transfinite Line {2, 4} = 60 Using Progression 1;

Transfinite Line {1,3} = 80 Using Bump 0.1;
Transfinite Line {2, 4} = 240 Using Progression 1;
Transfinite Surface {6};// Alternated;
*/

Physical Line("Gamma4") = {5,6}; //inlet
Physical Line("Gamma3") = {1}; //bottom
Physical Line("Gamma2") = {2,3}; //outlet
Physical Line("Gamma1") = {4}; //top
Physical Surface("Omega") = {9,11};

