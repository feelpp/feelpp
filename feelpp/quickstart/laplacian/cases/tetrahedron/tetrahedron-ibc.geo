h=0.1;
Point(1) = {-1, -1, -1, h};
Point(2) = {1, -1, -1, h};
Point(3) = {-1, 1, -1, h};
Point(4) = {-1, -1, 1, h};

Line(1) = {2, 3};
Line(2) = {3, 1};
Line(3) = {1, 2};

Line(4) = {1, 4};
Line(5) = {2, 4};
Line(6) = {3, 4};

//+                                                                                                                                                           
Line Loop(1) = {6, -5, 1};                                                                                                                                    
//+                                                                                                                                                           
Plane Surface(1) = {1};                                                                                                                                       
//+                                                                                                                                                           
Line Loop(2) = {2, 4, -6};                                                                                                                                    
//+                                                                                                                                                           
Plane Surface(2) = {2};                                                                                                                                       
//+                                                                                                                                                           
Line Loop(3) = {3, 5, -4};                                                                                                                                    
//+                                                                                                                                                           
Plane Surface(3) = {3};                                                                                                                                       
//+                                                                                                                                                           
Line Loop(4) = {3, 1, 2};                                                                                                                                     
//+                                                                                                                                                           
Plane Surface(4) = {4};                                                                                                                                       
//+                                                                                                                                                           
Surface Loop(1) = {2, 4, 3, 1};
Physical Surface("Dirichlet") = {1};
Physical Surface("Ibc") = {2};
Physical Surface("Neumann") = {3};
Physical Surface("Robin") = {4};
//+                                                                                                                                                           
Volume(0) = {1};                                                                                                                                              
Physical Volume("Tetrahedron") = {0};    



