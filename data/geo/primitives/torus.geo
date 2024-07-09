NL=20;
h = 0.005;
R=0.2;
C=2;

R1=C-R/2;
R2=C+R/2;
shape = 1; // square=1 circle=2
alphadiff=Pi/6;
alpha1=Pi/2;
alpha2=alpha1+alphadiff;
If ( shape == 1 )
  Point(1) = {0,C-R/2,-R/2,h};
  Point(2) = {0,C+R/2,-R/2,h};
  Point(3) = {0,C+R/2,+R/2,h};
  Point(4) = {0,C-R/2,+R/2,h};

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

Extrude {  {0,0,1}, {0,0,0}, alpha2-alpha1 }  {Surface{6}; Layers{C*(alpha2-alpha1)*R/h}; }  

Physical Surface("inlet") = {28};
Physical Surface("outlet") = {6};
Physical Surface("wall") = {19, 15, 27, 23};

Physical Volume("fluid") = {1};


