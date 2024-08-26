A = 30;
B = 10;
h=0.4;

y = 10;

Point(5) = {0, y, 0, h};
Point(6) = {0-A*0.5, y, 0, h};
Point(7) = {0+A*0.5, y, 0, h};
Point(8) = {0, y+0.5*B, 0, h};
Point(9) = {0, y-0.5*B, 0, h};

Ellipse(5) = {7, 5, 8};
Ellipse(6) = {8, 5, 6};
Ellipse(7) = {6, 5, 9};
Ellipse(8) = {9, 5, 7};

Curve Loop(2) = {5, 6, 7, 8};
Surface(1) = {2};

Rotate {{0, 0, 1}, {0, y, 0}, Pi/8} {
  Surface{1}; 
}

Physical Curve("Wall", 4) = {5, 6, 7, 8};
Physical Surface("Caoutchouc", 6) = {1};

Point(100) = {-25, 0, 0, h};
Point(200) = {-25, -1, 0, h};
Point(300) = {25, -1, 0, h};
Point(400) = {25, 0, 0, h};

Line(100) = {300, 200};
Line(200) = {200, 100};
Line(300) = {100, 400};
Line(400) = {400, 300};

Curve Loop(100) = {400, 100, 200, 300};
Plane Surface(100) = {100};

Physical Curve("Obs1", 500) = {400, 100, 200, 300};
Physical Surface("obstacle", 600) = {100};//+