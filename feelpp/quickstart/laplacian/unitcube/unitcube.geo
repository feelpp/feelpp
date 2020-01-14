// Gmsh project created on Tue Jul 23 09:33:33 2019                                                                                        
SetFactory("OpenCASCADE");                                                                                                                 
//+                                                                                                                                        
Box(1) = {0, 0, 0, 1, 1, 1};                                                                                                               
//+                                                                                                                                        
Physical Surface("Dirichlet") = {1, 5, 2, 3, 4};                                                                                           
//+                                                                                                                                        
Physical Volume("Omega") = {1};                                                                                                            
                                                                                                                                           
       
