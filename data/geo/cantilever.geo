h = 0.015;

height = 1;
length = 1.5;
offset = 0.5;
tip = 0.2;
ratio = 0.8;
Point(1) = {0, 0, 0, h};
Point(2) = {0, height, 0, h};
Point(3) = {0, height-offset/2, 0, h};
Point(4) = {0, offset/2, 0, h};
Point(5) = {length, (height+tip)/2, 0, h};
Point(6) = {length, (height-tip)/2, 0, h};
Point(7) = {0, height/2, 0, h};
Circle(1) = {4, 7, 3};
Point(8) = {2.5, 0.1, -0, h};
Point(9) = {2.5, 0.9, -0, h};
Ellipse(3) = {2, 1, 8, 5};
Ellipse(4) = {1, 2, 9, 6};
Line(5) = {3, 2};
Line(6) = {4, 1};
Line(7) = {6, 5};
Dilate {{0.95, 0.525, 0}, 0.45} {
  Duplicata { Line{3}; }
}
Dilate {{0.95, 0.475, 0}, 0.45} {
  Duplicata { Line{4}; }
}
Translate {0, -0.1, 0} {
  Point{15};
}
Translate {0, 0.1, 0} {
  Point{11};
}
Line(10) = {13, 17};
Ellipse(11) = {14, 7, 11, 10};
Line Loop(12) = {3, -7, -4, -6, 1, 5};
Line Loop(13) = {8, 10, -9, 11};
Plane Surface(14) = {12, 13};
Physical Line("wall") = {5, 6};
Physical Line("tip") = {7};
Physical Surface("cantilever")={14};
