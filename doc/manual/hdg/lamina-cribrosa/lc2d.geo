// cl__1 = 1e+22;


h=0.005; 

Mesh.CharacteristicLengthFactor=h;
Mesh.CharacteristicLengthMax=0.1;

Point(1) = {0.095, 0, -0.015, 1e+22};
Point(2) = {0, 0, -0.015, 1e+22};
Point(3) = {0.04750000000000001, 0.08227241335952167, -0.015, 1e+22};
Point(4) = {-0.08227241335952169, -0.04749999999999997, -0.015, 1e+22};
Point(5) = {0.008, 0, -0.015, 1e+22};
Point(6) = {0, 0.008, -0.015, 1e+22};
Point(7) = {-0.008, 0, -0.015, 1e+22};
Point(8) = {0, -0.008, -0.015, 1e+22};
Circle(1) = {1, 2, 3};
Circle(2) = {3, 2, 4};
Circle(3) = {4, 2, 1};
Circle(5) = {5, 2, 6};
Circle(6) = {6, 2, 7};
Circle(7) = {7, 2, 8};
Circle(8) = {8, 2, 5};
Line Loop(11) = {2, 3, 1, -7, -6, -5, -8};
Plane Surface(11) = {11};

Physical Line("top") = {5, 6, 7, 8}; // Dirichlet
Physical Line("bottom") = {1, 2, 3}; // Integral
Physical Surface("omega") = {11}; 
