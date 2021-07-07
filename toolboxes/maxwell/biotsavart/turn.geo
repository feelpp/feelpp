// Define Main params
h = 0.1;
r1=1;
r2=2;
he=2.5;
eps=0.1;
theta1=Asin( eps/(2*r1) );
theta2=Asin( eps/(2*r2) );

// 1st quarter
Point(1) = {0, 0, -he, h};

Point(2) = {r1*Cos(theta1), eps/2., -he, h};
Point(3) = {r2*Cos(theta2), eps/2., -he, h};
Point(4) = {0, r1, -he, h};
Point(5) = {0, r2, -he, h};
Point(6) = {-r1, 0, -he, h};
Point(7) = {-r2, 0, -he, h};
Point(8) = {0, -r1, -he, h};
Point(9) = {0, -r2, -he, h};
Point(10) = {r1*Cos(-theta1), -eps/2., -he, h};
Point(11) = {r2*Cos(-theta2), -eps/2., -he, h};

Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 6};
Circle(3) = {6, 1, 8};
Circle(4) = {8, 1, 10};

Circle(5) = {3, 1, 5};
Circle(6) = {5, 1, 7};
Circle(7) = {7, 1, 9};
Circle(8) = {9, 1, 11};

Line(9) = {2, 3};
Line(10) = {10, 11};

dL=newl; Line Loop(dL) = {1:4, 10, -8, -7, -6, -5, -9};
S=news; Plane Surface(S) = {dL};

out[] = Extrude {0,0,2*he} {Surface{S};};

Physical Volume("omega") = {out[1]};
Physical Surface("top") = {out[0]};
Physical Surface("bottom") = {S};
Physical Surface("Rint") = {out[2], out[3], out[4], out[5]};
Physical Surface("Rext") = {out[7], out[8], out[9], out[10]};
Physical Surface("V0") = {out[6]};
Physical Surface("V1") = {out[11]};
