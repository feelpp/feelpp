h=0.03;
Hx=1;
Hy=1;
Point(1) = {0., 0., 0., h};
Point(2) = {0, Hy, 0., h};
Point(3) = {Hx, Hy, 0., h};
Point(4) = {Hx, 0, 0., h};
Point(5) = {Hx/2, Hy/2, 0., h};
Point(6) = {0, Hy/2, 0., h};
Point(7) = {Hx, Hy/2, 0., h};
Point(8) = {Hx/2, 0, 0., h};
Point(9) = {Hx/2, Hy, 0., h};

Line(1) = {1, 8};
Line(2) = {8, 4};
Line(3) = {4, 7};
Line(4) = {7, 3};
Line(5) = {3, 9};
Line(6) = {9, 2};
Line(7) = {2, 6};
Line(8) = {6, 1};
Line(9) = {6, 5};
Line(10) = {5, 8};
Line(11) = {5, 7};
Line(12) = {5, 9};
Curve Loop(1) = {6, 7, 9, 12};
Plane Surface(1) = {1};
Curve Loop(2) = {5, -12, 11, 4};
Plane Surface(2) = {2};
Curve Loop(3) = {9, 10, -1, -8};
Plane Surface(3) = {3};
Curve Loop(4) = {11, -3, -2, -10};
Plane Surface(4) = {4};


Physical Surface("mat1") = {1};
Physical Surface("mat2") = {2};
Physical Surface("mat3") = {3};
Physical Surface("mat4") = {4};

Physical Surface("mat1_2") = {1,2};
Physical Surface("mat1_2_3") = {1,2,3};
Physical Surface("mat1_2_3_4") = {1,2,3,4};


Physical Line("gamma_x0_mat3") = {8};
Physical Line("gamma_x0_mat1") = {7};
Physical Line("gamma_xH_mat2") = {4};
Physical Line("gamma_xH_mat4") = {3};

Physical Line("gamma_y0_mat3") = {1};
Physical Line("gamma_y0_mat4") = {2};
Physical Line("gamma_yH_mat1") = {6};
Physical Line("gamma_yH_mat2") = {5};

Physical Line("interface_mat1_mat2") = {12};
Physical Line("interface_mat1_mat3") = {9};
Physical Line("interface_mat3_mat4") = {10};
Physical Line("interface_mat2_mat4") = {11};

Physical Line("gamma_x0_mat1_3") = {7,8};
Physical Line("gamma_all") = {1,2,3,4,5,6,7,8};
Physical Line("interface_all") = {9,10,11,12};

Physical Point("point_corner_00") = {1};
Physical Point("point_corner_all") = {1,2,3,4};
Physical Point("point_geo_all") = {1,2,3,4,5,6,7,8,9};
