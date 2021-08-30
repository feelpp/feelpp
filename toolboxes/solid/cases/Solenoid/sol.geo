h=0.1;
r_int=1;
r_ext=2;
z=20;
Point(1) = {0, 0, 0, h};                               
Point(2) = {r_int, 0, 0, h};                           
Point(3) = {-r_int, 0, 0, h};                          
Point(4) = {0, r_int, 0, h};                           
Point(5) = {0, -r_int, 0, h};                          
Point(7) = {r_ext, 0, 0, h};                           
Point(8) = {-r_ext, 0, 0, h};                          
Point(9) = {0, r_ext, 0, h};                           
Point(10) = {0, -r_ext, 0, h};                         
Circle(1) = {2, 1, 4};                                  
Circle(2) = {4, 1, 3};                                  
Circle(3) = {3, 1, 5};                                  
Circle(4) = {5, 1, 2};                                  
Circle(5) = {7, 1, 9};                                  
Circle(6) = {9, 1, 8};                                  
Circle(7) = {8, 1, 10};                                 
Circle(8) = {10, 1, 7};                                 
Line Loop(9) = {5, 6, 7, 8};                            
Line Loop(10) = {1, 2, 3, 4};                           
Plane Surface(11) = {9, 10};                            
Extrude {0, 0, z/2} {                                   
  Surface{11};                                          
}                                                       
Extrude {0, 0, z/2} {                                   
  Surface{53};                                          
}                                                       
Physical Surface("Bottom") = {11};                    
Physical Surface("Middle") = {53};                    
Physical Surface("Top") = {95};                       
Physical Surface("Rext") = {24,28,32,36,66,70,74,78}; 
Physical Surface("Rint") = {40,44,48,52,82,86,90,94}; 
Physical Volume("Omega") = {1,2};
