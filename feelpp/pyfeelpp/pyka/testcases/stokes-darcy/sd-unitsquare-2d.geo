h=0.1;
//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {0.5, 0, 0, h};
//+
Point(3) = {1, 0, 0, h};
//+
Point(4) = {1, 1, 0, h};
//+
Point(5) = {0.5, 1, 0, h};
//+
Point(6) = {0, 1, 0, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Line(7) = {2, 5};
//+
Line Loop(1) = {6, 1, 7, 5};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {7, -4, -3, -2};
//+
Plane Surface(2) = {2};
//+
Physical Line("Gamma") = {7};
//+
Physical Line("Gammaf-bottom") = {1};
//+
Physical Line("Gammaf-top") = {5};
//+
Physical Line("Gammaf-side") = {6};
//+
Physical Line("Gammap-side") = {3};
//+
Physical Line("Gammap-top") = {4};
//+
Physical Line("Gammap-bottom") = {2};
//+
Physical Surface("Omegap") = {2};
//+
Physical Surface("Omegaf") = {1};
