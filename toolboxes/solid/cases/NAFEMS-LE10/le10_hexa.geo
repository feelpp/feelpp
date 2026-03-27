h = 0.1;

Point(1) = {0,0,0,h};
Point(2) = {0,1,0,h};
Point(3) = {2,0,0,h};
Point(4) = {3.25,0,0,h};
Point(5) = {0,2.75,0,h};

Ellipse(1) = {2,1,2,3}; 
Line(2) = {3,4};
Ellipse(3) = {4,1,4,5};
Line(4) = {5,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(1) = {5};

Transfinite Curve {1, 3} = 25*2; // Bords courbes (arcs)
Transfinite Curve {2, 4} = 15*2; // Bords droits (rayons)

Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, 0, 0.3} {
  Surface{1}; Layers{3}; Recombine;
}
Extrude {0, 0, 0.3} {
  Surface{27}; Layers{3}; Recombine;
}

Physical Line( "EE" ) = {9};
Physical Surface( "Bottom" ) = {1};                                                                                                                                                                                                                                                                                                                                                                                       
Physical Surface( "Top" ) = {49};
Physical Surface( "DCDC" ) = {18,40};
Physical Surface( "ABAB" ) = {26,48};
Physical Surface( "BCBC" ) = {22,44};
Physical Surface( "ADAD" ) = {14,36};
Physical Volume( "Omega" ) = {1,2};