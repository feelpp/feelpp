h=0.01;
H=0.2;
L=4*H;

Point (1) = {0,0,0,h} ;
Point (2) = {L,0,0,h} ;
Point (3) = {L,H,0,h} ;
Point (4) = {0,H,0,h} ;

Point (5) = {L/2,0,0,h/2};
Point (6) = {L/2,H/2,0,h/2} ;
Point (7) = {L/2,H,0,h/2} ;
Point (8) = {0,H/2,0,h} ;
Point (9) = {L,H/2,0,h} ;


Line (1)={1,5};
Line (2)={5,6};
Line (3)={6,8};
Line (4)={8,1};

Line (5) = {5,6};
Line (6) = {6,9};
Line (7) = {9,6};
Line (8) = {6,5};

Line (9) = {9,3};
Line (10) = {3,7};
Line (11) = {7,6};
Line (12) = {7,4};
Line (13) = {4,8};
Line (14) = {5,2};
Line (15) = {2,9};

Curve Loop (1) = {1,2,3,4} ;
Plane Surface (1) = {1} ;

Curve Loop (2) = {14, 15, -6, -2};
Plane Surface (2) = {2};

Curve Loop (3) = {6,9,10,11};
Plane Surface (3) = {3};

Curve Loop (4) = {-3,-11,12,13};
Plane Surface (4) = {4};

Physical Surface("Omega1") = {1};
Physical Surface("Omega2") = {2};
Physical Surface("Omega3") = {3};
Physical Surface("Omega4") = {4};
