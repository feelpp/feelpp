h=0.03;

x1=0.0;
x2=0.7;
x3=0.8;
x4=1.0;

y1=0.0;
y2=0.2;
y3=0.6;
y4=0.7;
y5=0.8;
y6=1.0;

Point(1)={x1,y1,0,h};
Point(2)={x4,y1,0,h};
Point(3)={x4,y6,0,h};
Point(4)={x1,y6,0,h};
Point(5)={x1,y5,0,h};
Point(6)={x1,y4,0,h};

Point( 7)={x2,y2,0,h};
Point( 8)={x3,y2,0,h};
Point( 9)={x3,y3,0,h};
Point(10)={x2,y3,0,h};

Line( 1) = {1,2};
Line( 2) = {2,3};
Line( 3) = {3,4};
Line( 4) = {4,5};
Line( 5) = {5,6};
Line( 6) = {6,1};

Line( 7) = { 7, 8};
Line( 8) = { 8, 9};
Line( 9) = { 9,10};
Line(10) = {10, 7};

Line Loop(1)={1,2,3,4,5,6};
Line Loop(2)={7,8,9,10};
Line Loop(3)={5};

Plane Surface(1)={1,-2};
Plane Surface(2)={2};

Physical Line("border")={6,1,2,3,4};
Physical Line("border_steak")={7,8,9,10};
Physical Line("source")={5};

Physical Surface("oven")={1};
Physical Surface("steak")={2};
