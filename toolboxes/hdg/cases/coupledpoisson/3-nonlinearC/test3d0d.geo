h=0.1;
Point(1) = {0, 0, 0, h};
Point(2) = {0, 1, 1, h};
Point(3) = {0, 0, 1, h};
Point(4) = {0, 1, 0, h};
//+
Line(1) = {4, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 1};
//+
Line(4) = {1, 4};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Extrude {2, 0, 0} {
  Surface{1}; 
}
//+
Physical Volume("omega") = {1};
//+
Physical Surface("bottom") = {1};
//+
Physical Surface("top") = {26};
//+
Physical Surface("lateral") = {17, 21, 13, 25};
