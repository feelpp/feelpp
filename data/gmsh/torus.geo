h = 0.05;
R=0.2;
tier=1; // one third of a torus
shape=2; // square=1 circle=2
If ( shape == 1 )
  Point(1) = {-R,-R,0,h};
  Point(2) = {R,-R,0,h};
  Point(3) = {R,R,0,h};
  Point(4) = {-R,R,0,h};
  Line(1) = {1,2};
  Line(2) = {2,3};
  Line(3) = {3,4};
  Line(4) = {4,1};
EndIf
If ( shape == 2 )
  Point(1) = {0,0,0,h};
  Point(2) = {0,R,0,h};
  Point(3) = {0,0,R,h};
  Point(4) = {0,0,-R,h};
  Point(5) = {0,-R,0,h};

  Circle(1) = {2,1,3};
  Circle(2) = {3,1,5};
  Circle(3) = {5,1,4};
  Circle(4) = {4,1,2};
EndIf

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
If ( tier == 1 )
  Extrude Surface {6, {0,0,1}, {0,-2,0}, 2*Pi/3}{Layers{1./h};};
  
  Physical Surface("inlet") = {28};
  Physical Surface("outlet") = {6};
  Physical Surface("wall") = {19, 15, 27, 23};
  
  Physical Volume("fluid") = {1};
EndIf
If ( tier == 2 )
  Extrude Surface {28, {0,0,1}, {0,-2,0}, 2*Pi/3}{Layers{1/.h};};
EndIf
If ( tier == 3 )
  Extrude Surface {50, {0,0,1}, {0,-2,0}, 2*Pi/3}{Layers{1/.h};};
EndIf




