/** -*- mode: gmsh -*-
  This is the step defined in 

  Noundary conditions in approximate commutator preconditioners for the navier-stokes equations 
  Howard C. Elman and Ray S. Tuminaro
  2009

  Author : Christophe Prud'homme      
*/
dim=3;
h=0.1;
L=10;
W=1;
Point(1) = {-1, -1, 0, h};
Point(2) = {0, -1, 0, h};
Point(3) = {L, -1, 0, h};
Point(4) = {L, 0, 0, h};
Point(5) = {L, 1, 0, h};
Point(6) = {0, 1, 0, h};
Point(7) = {-1, 1, 0, h};
Point(8) = {-1, 0, 0, h};
Point(9) = {0, 0, 0, h};
Line(1) = {8, 9};
Line(2) = {9, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line Loop(9) = {6, 7, 8, 1, 2, 3, 4, 5};
Plane Surface(10) = {9};

If ( dim == 2 )
  Physical Line("inlet") = {8};
  Physical Line("wall") = {6, 7, 1, 2, 3};
  Physical Line("outlet") = {5, 4};
  Physical Surface("fluid") = {10};
EndIf 
If (dim == 3 )
  Extrude {0, 0, W} {
    Surface{10};
  }
  Surface Loop(53) = {52, 23, 10, 27, 31, 35, 39, 43, 47, 51};
  Volume(54) = {53};
  Physical Surface("inlet") = {31};
  Physical Surface("outlet") = {51, 47};
  Physical Surface("wall") = {23, 10, 52, 43, 27, 39, 35};
  Physical Volume("fluid") = {1};
EndIf
