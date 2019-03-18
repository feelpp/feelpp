h=0.1;
//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {1, 0, 0, h};
//+
Point(3) = {1, 1, 0, h};
//+
Point(4) = {0, 1, 0, h};
//+
Point(5) = {0, 0, 0.5, h};
//+
Point(6) = {1, 0, 0.5, h};
//+
Point(7) = {1, 1, 0.5, h};
//+
Point(8) = {0, 1, 0.5, h};
//+
Point(9) = {0, 0, 1, h};
//+
Point(10) = {1, 0, 1, h};
//+
Point(11) = {1, 1, 1, h};
//+
Point(12) = {0, 1, 1, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 5};
//+
Line(8) = {5, 6};
//+
Line(9) = {10, 11};
//+
Line(10) = {11, 12};
//+
Line(11) = {12, 9};
//+
Line(12) = {9, 10};
//+
Line(13) = {2, 6};
//+
Line(14) = {6, 10};
//+
Line(15) = {3, 7};
//+
Line(16) = {7, 11};
//+
Line(17) = {4, 8};
//+
Line(18) = {8, 12};
//+
Line(19) = {1, 5};
//+
Line(20) = {5, 9};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {7, 8, 5, 6};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {12, 9, 10, 11};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {17, -6, -15, 3};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {6, 18, -10, -16};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {2, 15, -5, -13};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {5, 16, -9, -14};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {12, -14, -8, 20};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {8, -13, -1, 19};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {4, 19, -7, -17};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {7, 20, -11, -18};
//+
Plane Surface(11) = {11};
//+
Surface Loop(1) = {1, 10, 9, 6, 4, 2};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {3, 8, 7, 5, 11, 2};
//+
Volume(2) = {2};
//+
Physical Surface("Gamma") = {2};
//+
Physical Surface("Gammap-bottom") = {1};
//+
Physical Surface("Gammap-side-south") = {9};
//+
Physical Surface("Gammap-side-east") = {6};
//+
Physical Surface("Gammap-side-north") = {4};
//+
Physical Surface("Gammap-side-west") = {10};
//+
Physical Surface("Gammaf-side-west") = {11};
//+
Physical Surface("Gammaf-side-south") = {8};
//+
Physical Surface("Gammaf-side-east") = {7};
//+
Physical Surface("Gammaf-side-north") = {5};
//+
Physical Surface("Gammaf-top") = {3};
//+
Physical Volume("Gammap") = {1};
//+
Physical Volume("Gammaf") = {2};
