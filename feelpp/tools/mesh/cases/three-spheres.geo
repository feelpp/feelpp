// Mesh where the three spheres are holes in the fluid; they are not filled
RSphere = 1.;
lcSphere = .25;

RDom = 10;
lcDom = 4.;


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

Right = 10*RSphere;

Point(101) = {Right,0,0,lcSphere};
Point(102) = {RSphere+Right,0,0,lcSphere};
Point(103) = {Right,RSphere,0,lcSphere};
Circle(101) = {102,101,103};
Point(104) = {-RSphere+ Right,0,0.0,lcSphere};
Point(105) = {Right,-RSphere,0.0,lcSphere};
Circle(102) = {103,101,104};
Circle(103) = {104,101,105};
Circle(104) = {105,101,102};
Point(106) = {Right,0,-RSphere,lcSphere};
Point(107) = {Right,0,RSphere,lcSphere};
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
Volume(2222) = {129};

Left = -10*RSphere;

Point(201) = {Left,0,0,lcSphere};
Point(202) = {RSphere+Left,0,0,lcSphere};
Point(203) = {Left,RSphere,0,lcSphere};
Circle(201) = {202,201,203};
Point(204) = {-RSphere+ Left,0,0.0,lcSphere};
Point(205) = {Left,-RSphere,0.0,lcSphere};
Circle(202) = {203,201,204};
Circle(203) = {204,201,205};
Circle(204) = {205,201,202};
Point(206) = {Left,0,-RSphere,lcSphere};
Point(207) = {Left,0,RSphere,lcSphere};
Circle(205) = {203,201,206};
Circle(206) = {206,201,205};
Circle(207) = {205,201,207};
Circle(208) = {207,201,203};
Circle(209) = {202,201,207};
Circle(210) = {207,201,204};
Circle(211) = {204,201,206};
Circle(212) = {206,201,202};
Line Loop(213) = {202,208,-210};
Surface(214) = {213};
Line Loop(215) = {210,203,207};
Surface(216) = {215};
Line Loop(217) = {-208,-209,201};
Surface(218) = {217};
Line Loop(219) = {-211,-202,205};
Surface(220) = {219};
Line Loop(221) = {-205,-212,-201};
Surface(222) = {221};
Line Loop(223) = {-203,211,206};
Surface(224) = {223};
Line Loop(225) = {-207,204,209};
Surface(226) = {225};
Line Loop(227) = {-204,212,-206};
Surface(228) = {227};
Surface Loop(229) = {228,226,216,214,220,224,222,218};
Volume(3333) = {229};

HalfSide = 5*Right;
HalfHeight = 2*Right;

Point(401) = {HalfSide, HalfHeight, HalfHeight,lcDom};
Point(402) = {HalfSide, HalfHeight,-HalfHeight,lcDom};
Point(403) = {HalfSide,-HalfHeight, HalfHeight,lcDom};
Point(404) = {HalfSide,-HalfHeight,-HalfHeight,lcDom};
Point(405) = {-HalfSide, -HalfHeight, -HalfHeight,lcDom};
Point(406) = {-HalfSide,HalfHeight, -HalfHeight,lcDom};
Point(407) = {-HalfSide,-HalfHeight, HalfHeight,lcDom};
Point(408) = {-HalfSide, HalfHeight, HalfHeight,lcDom};


//+
Line(413) = {408, 406};
//+
Line(414) = {406, 405};
//+
Line(415) = {405, 407};
//+
Line(416) = {407, 408};
//+
Line(417) = {408, 401};
//+
Line(418) = {401, 403};
//+
Line(419) = {403, 404};
//+
Line(420) = {404, 405};
//+
Line(421) = {403, 407};
//+
Line(422) = {404, 402};
//+
Line(423) = {402, 401};
//+
Line(424) = {406, 402};



//+
Curve Loop(428) = {416, 413, 414, 415};
//+
Plane Surface(429) = {428};
//+
Curve Loop(429) = {418, 419, 422, 423};
//+
Plane Surface(430) = {429};
//+
Curve Loop(430) = {421, -415, -420, -419};
//+
Plane Surface(431) = {430};
//+
Curve Loop(431) = {417, -423, -424, -413};
//+
Plane Surface(432) = {431};
//+
Curve Loop(432) = {416, 417, 418, 421};
//+
Plane Surface(433) = {432};
//+
Curve Loop(433) = {420, -414, 424, -422};
//+
Plane Surface(434) = {433};

Physical Surface("CenterSphere") = {28,26,16,14,20,24,22,18};
Physical Surface("RightSphere") = {128,126,116,114,120,124,122,118};
Physical Surface("LeftSphere") = {228,226,216,214,220,224,222,218};
Physical Surface("Wall") = {429, 431, 433, 432, 434, 430};
//+
//Physical Volume("Body") = {1111,2222,3333};
Surface Loop(230) = {431, 433, 429, 432, 430, 434};
//+
Volume(4444) = {29, 129, 229, 230};
//+
Physical Volume("Fluid") = {4444};