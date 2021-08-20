h = 0.1;
nNodesByLine = 20;
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
Point(1) = {xmin,ymin,0.0,h};
Point(2) = {xmax,ymin+0,0.0,h};
Point(3) = {xmax+0,ymax,0.0,h};
Point(4) = {xmin+0,ymax+0,0.0,h};
Line(1) = {4,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Physical Point("Top-Left") = {4};
Physical Point("Top-Right") = {3};
Physical Point("Bottom-Right") = {2};
Physical Point("Bottom-Left") = {1};
Physical Line("Left") = {1};
Physical Line("Bottom") = {2};
Physical Line("Right") = {3};
Physical Line("Top") = {4};
Physical Surface("Omega") = {6};

Transfinite Line {1,2,3,4} = nNodesByLine Using Progression 1;
Transfinite Surface {6} Right;//Left;//Right; //  Left | Right | Alternate | AlternateRight | AlternateLeft >
