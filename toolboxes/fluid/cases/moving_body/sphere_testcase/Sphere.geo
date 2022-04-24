
RSphere = 0.1;
lcSphere =0.0125;
RDom = 10;
lcDom = 1.5;


Point(1) = {0,0,0,lcSphere};
Point(2) = {RSphere,0,0,lcSphere};
Point(3) = {0,RSphere,0,lcSphere};
Circle(1) = {2,1,3};
Point(4) = {-RSphere,0,0.0,lcSphere};
Point(5) = {0,-RSphere,0.0,lcSphere};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Point(6) = {0,0,-RSphere,lcSphere};
Point(7) = {0,0,RSphere,lcSphere};
Circle(5) = {3,1,6};
Circle(6) = {6,1,5};
Circle(7) = {5,1,7};
Circle(8) = {7,1,3};
Circle(9) = {2,1,7};
Circle(10) = {7,1,4};
Circle(11) = {4,1,6};
Circle(12) = {6,1,2};
Line Loop(13) = {2,8,-10};
Surface(14) = {13};
Line Loop(15) = {10,3,7};
Surface(16) = {15};
Line Loop(17) = {-8,-9,1};
Surface(18) = {17};
Line Loop(19) = {-11,-2,5};
Surface(20) = {19};
Line Loop(21) = {-5,-12,-1};
Surface(22) = {21};
Line Loop(23) = {-3,11,6};
Surface(24) = {23};
Line Loop(25) = {-7,4,9};
Surface(26) = {25};
Line Loop(27) = {-4,12,-6};
Surface(28) = {27};
Surface Loop(29) = {28,26,16,14,20,24,22,18};
Volume(1111) = {29};

Point(101) = {0,0,0,lcDom};
Point(102) = {RDom,0,0,lcDom};
Point(103) = {0,RDom,0,lcDom};
Circle(101) = {102,101,103};
Point(104) = {-RDom,0,0.0,lcDom};
Point(105) = {0,-RDom,0.0,lcDom};
Circle(102) = {103,101,104};
Circle(103) = {104,101,105};
Circle(104) = {105,101,102};
Point(106) = {0,0,-RDom,lcDom};
Point(107) = {0,0,RDom,lcDom};
Circle(105) = {103,101,106};
Circle(106) = {106,101,105};
Circle(107) = {105,101,107};
Circle(108) = {107,101,103};
Circle(109) = {102,101,107};
Circle(110) = {107,101,104};
Circle(111) = {104,101,106};
Circle(112) = {106,101,102};
Line Loop(113) = {102,108,-110};
Surface(114) = {113};
Line Loop(115) = {110,103,107};
Surface(116) = {115};
Line Loop(117) = {-108,-109,101};
Surface(118) = {117};
Line Loop(119) = {-111,-102,105};
Surface(120) = {119};
Line Loop(121) = {-105,-112,-101};
Surface(122) = {121};
Line Loop(123) = {-103,111,106};
Surface(124) = {123};
Line Loop(125) = {-107,104,109};
Surface(126) = {125};
Line Loop(127) = {-104,112,-106};
Surface(128) = {127};
Surface Loop(129) = {128,126,116,114,120,124,122,118};

Volume(200) = {129,29};


Physical Surface("Sphere") = {28,26,16,14,20,24,22,18};
Physical Surface("Walls") = {114,118,120,122};
Physical Surface("Outlet") = {128,126,116,124};
Physical Volume("Fluid") = {200};
Physical Volume("SphereVol") = {1111};




