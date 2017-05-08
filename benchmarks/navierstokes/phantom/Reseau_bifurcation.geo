// Gmsh project created on Thu May 22 14:13:10 2014
h=200;
hmax = 0.4;
hmean=0.2;
hmin=0.1;
//Input begin
Point(1) = {0, 0, 0, hmax};
Point(2) = {0, 2.5, 0, hmax};
Point(3) = {2.5, 0, 0, hmax};
Point(4) = {-2.5, 0, 0, hmax};
Point(5) = {0, -2.5, 0, hmax};
Circle(1) = {5, 1, 3};
Circle(2) = {3, 1, 2};
Circle(3) = {2, 1, 4};
Circle(4) = {4, 1, 5};
//Input end
Point(6) = {0, 2.5, 40, hmax};
Point(7) = {2.5, 0, 40, hmax};
Point(8) = {-2.5, 0, 40, hmax};
Point(9) = {0, -2.5, 40, hmax};
Point(10) = {0, 0, 40, hmean};
Circle(9) = {9, 10, 7};
Circle(10) = {7, 10, 6};
Circle(11) = {6, 10, 8};
Circle(12) = {8, 10, 9};
//Little input
Point(11) = {1.5, 0, 40, hmean};
Point(12) = {-1.5, 0, 40, hmean};
Point(13) = {0, 1.5, 40, hmean};
Point(14) = {0, -1.5, 40, hmean};
Circle(14) = {11, 10, 13};
Circle(15) = {13, 10, 12};
Circle(16) = {12, 10, 14};
Circle(17) = {14, 10, 11};
//First branch
Point(15) = {0, 6.61, 69.88, hmean};
Point(16) = {0, 3.61, 69.88, hmean};
Point(17) = {0, 5.11, 69.88, hmean};
Point(18) = {1.5, 5.11, 69.88, hmean};
Point(19) = {-1.5, 5.11, 69.88, hmean};
Circle(18) = {18, 17, 15};
Circle(19) = {15, 17, 19};
Circle(20) = {19, 17, 16};
Circle(21) = {16, 17, 18};
//Second branch
Point(20) = {0, -6.11, 69.88, hmin};
Point(21) = {0, -4.11, 69.88, hmin};
Point(22) = {0, -5.11, 69.88, hmin};
Point(23) = {-1, -5.11, 69.88, hmin};
Point(24) = {1, -5.11, 69.88, hmin};
Circle(25) = {24, 22, 21};
Circle(26) = {21, 22, 23};
Circle(27) = {23, 22, 20};
Circle(28) = {20, 22, 24};
//Little output
Point(25) = {1.5, 0, 99.77, hmean};
Point(26) = {-1.5, 0, 99.77, hmean};
Point(27) = {0, 1.5, 99.77, hmean};
Point(28) = {0, -1.5, 99.77, hmean};
Point(29) = {0, 0, 99.77, hmean};
Circle(29) = {25, 29, 27};
Circle(30) = {27, 29, 26};
Circle(31) = {26, 29, 28};
Circle(32) = {28, 29, 25};
//Output begin
Point(30) = {2, 0, 99.77, hmax};
Point(31) = {-2, 0, 99.77, hmax};
Point(32) = {0, 2, 99.77, hmax};
Point(33) = {0, -2, 99.77, hmax};
Circle(33) = {30, 29, 32};
Circle(34) = {32, 29, 31};
Circle(35) = {31, 29, 33};
Circle(36) = {33, 29, 30};
//Output end
Point(34) = {2, 0, 139.77, hmax};
Point(35) = {-2, 0, 139.77, hmax};
Point(36) = {0, 2, 139.77, hmax};
Point(37) = {0, -2, 139.77, hmax};
Point(38) = {0, 0, 139.77, hmax};
Circle(37) = {34, 38, 36};
Circle(38) = {36, 38, 35};
Circle(39) = {35, 38, 37};
Circle(40) = {37, 38, 34};
//Easy lines
//Input
Line(41) = {2, 6};
Line(42) = {7, 3};
Line(43) = {4, 8};
Line(44) = {5, 9};
//Output
Line(45) = {32, 36};
Line(46) = {31, 35};
Line(47) = {30, 34};
Line(48) = {37, 33};
//First branch
//input
Line(49) = {13, 15};
Line(50) = {12, 19};
Line(51) = {11, 18};
Line(52) = {14, 16};
//output
Line(53) = {16, 28};
Line(54) = {15, 27};
Line(55) = {18, 25};
Line(56) = {19, 26};
//For second branch, we have to replot little circles
Point(39) = {1, 0, 40, h};
Point(40) = {-1, 0, 40, h};
Point(41) = {0, 1, 40, h};
Point(42) = {0, -1, 40, h};
Circle(57) = {39, 10, 41};
Circle(58) = {41, 10, 40};
Circle(59) = {40, 10, 42};
Circle(60) = {42, 10, 39};
Point(43) = {1, 0, 99.77, h};
Point(44) = {-1, 0, 99.77, h};
Point(45) = {0, 1, 99.77, h};
Point(46) = {0, -1, 99.77, h};
Circle(65) = {43, 29, 45};
Circle(66) = {45, 29, 44};
Circle(67) = {44, 29, 46};
Circle(68) = {46, 29, 43};
//Second branch
//input
Line(61) = {41, 21};
Line(62) = {42, 20};
Line(63) = {40, 23};
Line(64) = {39, 24};
//ouptut
Line(69) = {46, 20};
Line(70) = {43, 24};
Line(71) = {45, 21};
Line(72) = {23, 44};
//New points for intersection
//input
Point(47) = {0, -1.25, 41.4618395303, hmin};
Point(48) = {0, -0.25, 47.3091976517, hmin};
Delete {
  Line{52};
}
Line(73) = {16, 48};
Line(74) = {47, 14};
Delete {
  Line{62};
}
Delete {
  Line{61};
}
Line(75) = {21, 48};
Line(76) = {20, 47};
//output
Point(49) = {0, -1.25, 98.3076712329, hmin};
Point(50) = {0, -0.25, 92.4583561644, hmin};
Delete {
  Line{53};
}
Line(77) = {16, 50};
Line(78) = {49, 28};
Delete {
  Line{71};
}
Delete {
  Line{69};
}
Line(113) = {21, 50};
Line(114) = {20, 49};
//hard intersection of cylinders
Point(51) = {1, -1.11803398875, 40, hmean};
Point(52) = {1, 3.991966011, 69.88, hmean};
Point(53) = {-1, -1.11803398875, 40, hmean};
Point(54) = {-1, 3.991966011, 69.88, hmean};
Line(79) = {51, 52};
Line(80) = {53, 54};
Point(55) = {1, -1.11803398875, 99.77, hmean};
Point(56) = {-1, -1.11803398875, 99.77, hmean};
Line(81) = {55, 52};
Line(82) = {54, 56};
//input
Point(57) = {-1, -0.559016994389, 43.2687725621, hmean};
Point(58) = {+1, -0.559016994389, 43.2687725621, hmin};

Delete {
  Line{64};
}
Line(119) = {58, 24};
Delete {
  Line{63};
}
Line(83) = {57, 23};
//output
Point(59) = {-1, -0.559016994389, 96.5001334712, hmin};
Point(60) = {+1, -0.559016994389, 96.5001334712, hmin};
Delete {
  Line{72};
}
Line(84) = {23, 59};
Delete {
  Line{70};
}
Line(85) = {24, 60};
//spline second branch
//Spline(120) = {48, 58, 47, 57, 48};
//Spline(121) = {49, 60, 50, 59, 49};
//Delete some construction elements
Delete {
  Line{60};
}
Delete {
  Line{57, 58};
}
Delete {
  Line{59};
}
Delete {
  Point{10, 39};
}
Delete {
  Point{41};
}
Delete {
  Point{40};
}
Delete {
  Point{42};
}
Delete {
  Line{68};
}
Delete {
  Line{65};
}
Delete {
  Line{66};
}
Delete {
  Line{67};
}
Delete {
  Point{43};
}
Delete {
  Point{45};
}
Delete {
  Point{44};
}
Delete {
  Point{46};
}
//Split some lines
Delete {
  Line{21};
}
Circle(165) = {18, 17, 52};
Circle(166) = {52, 17, 16};
Delete {
  Line{20};
}
Circle(167) = {19, 17, 54};
Circle(168) = {54, 17, 16};
Delete {
  Line{17};
}
Circle(169) = {11, 10, 51};
Circle(170) = {51, 10, 14};
Delete {
  Line{16};
}
Circle(171) = {12, 10, 53};
Circle(172) = {53, 10, 14};
Delete {
  Line{32, 31};
}
Circle(179) = {26, 29, 56};
Circle(180) = {56, 29, 28};
Circle(181) = {28, 29, 55};
Circle(182) = {55, 29, 25};
//Split Line(121) {59, 49, 60, 50};
//Split Line(120) {48, 47, 58, 57};

Delete {
	Line{79, 80, 81, 82};
}
Line(191) = {54, 57};
Line(192) = {57, 53};
Line(193) = {52, 58};
Line(194) = {58, 51};
Line(195) = {52, 60};
Line(196) = {60, 55};
Line(197) = {54, 59};
Line(198) = {59, 56};

//Add circles
Point(61) = {0, 0.559016994361, 43.2687725621, hmean};
Point(62) = {1.5, 0.559016994361, 43.2687725621, hmean};
Point(63) = {-1.5, 0.559016994361, 43.2687725621, hmean};
Point(64) = {0, 2.05901699436, 43.2687725621, hmean};
Circle(199) = {58, 61, 62};
Circle(200) = {57, 61, 63};
Circle(201) = {62, 61, 64};
Circle(202) = {64, 61, 63};

Point(65) = {0, 0.559016994361, 96.5012274379, hmean};
Point(66) = {1.5, 0.559016994361, 96.5012274379, hmean};
Point(67) = {-1.5, 0.559016994361, 96.5012274379, hmean};
Point(68) = {0, 2.05901699436, 96.5012274379, hmean};
Circle(203) = {59, 65, 67};
Circle(204) = {67, 65, 68};
Circle(205) = {68, 65, 66};
Circle(206) = {66, 65, 60};

//Delete lines for split
Delete {
  Line{54};
}
Delete {
  Line{49};
}
Line(207) = {15, 68};
Line(208) = {68, 27};
Line(209) = {15, 64};
Line(210) = {64, 13};
Delete {
  Line{50};
}
Delete {
  Line{51};
}
Delete {
  Line{55};
}
Delete {
  Line{56};
}
Line(211) = {18, 62};
Line(212) = {62, 11};
Line(213) = {12, 63};
Line(214) = {63, 19};
Line(215) = {19, 67};
Line(216) = {67, 26};
Line(217) = {18, 66};
Line(218) = {66, 25};

//Improve branch connexion
Point(1005)={0.95918220939493,-0.65402249385266,42.846335298178,hmin};
Point(1006)={0.40284528306364,-0.29224258679214,46.493299852679,hmin};
Point(1007)={0.6395195417088,-0.34599342923419,45.592048818572,hmin};
Point(1008)={0.82756981496187,-0.41018410990994,44.679999836464,hmin};
Point(1009)={0.96648744559971,-0.49360432195332,43.760851081723,hmin};
Point(1011)={0.64184969510739,-0.95460949606404,42.026945386298,hmin};
Point(1013)={-0.64184969509349,-0.95460949607436,42.026945386276,hmin};
Point(1014)={-0.95918220939059,-0.65402249385891,42.846335298155,hmin};
Point(1015)={-0.96648744560175,-0.49360432195532,43.760851081705,hmin};
Point(1016)={-0.82756981496426,-0.41018410991096,44.679999836451,hmin};
Point(1017)={-0.63951954171064,-0.34599342923471,45.592048818564,hmin};
Point(1018)={-0.4028452830648,-0.29224258679235,46.493299852675,hmin};
Spline(187) = {58, 1005, 1011, 47};
Spline(188) = {47, 1013, 1014, 57};
Spline(189) = {57, 1015, 1016, 1017, 1018, 48};
Spline(190) = {48, 1006, 1007, 1008, 1009, 58};

Point(2001)={0.64191438925625,-0.95456145295052,97.742270537501,hmin};
Point(2002)={0.95920933492305,-0.6539834138538,96.922570632675,hmin};
Point(2003)={0.96647566310195,-0.49359273810579,96.007784492654,hmin};
Point(2004)={0.82755879213021,-0.41017939020716,95.088373245276,hmin};
Point(2005)={0.63951476697836,-0.34599210120464,94.176058991795,hmin};
Point(2006)={0.40284884641917,-0.29224323686685,93.274538744406,hmin};
Point(2007)={-0.40284884641059,-0.29224323686528,93.274538744378,hmin};
Point(2008)={-0.63951476697295,-0.34599210120314,94.176058991772,hmin};
Point(2009)={-0.82755879212678,-0.41017939020569,95.088373245258,hmin};
Point(2010)={-0.96647566310043,-0.49359273810429,96.00778449264,hmin};
Point(2011)={-0.95920933492461,-0.65398341385155,96.922570632667,hmin};
Point(2012)={-0.64191438925808,-0.95456145294916,97.742270537498,hmin};
Spline(183) = {60, 2003, 2004, 2005, 2006, 50};
Spline(184) = {50, 2007, 2008, 2009, 2010, 59};
Spline(185) = {59, 2011, 2012, 49};
Spline(186) = {49, 2001, 2002, 60};

//Surface
Line Loop(219) = {165, 193, 199, -211};
Ruled Surface(220) = {219};
Line Loop(221) = {18, 209, -201, -211};
Ruled Surface(222) = {221};
Line Loop(223) = {19, -214, -202, -209};
Ruled Surface(224) = {223};
Line Loop(225) = {214, 167, 191, 200};
Ruled Surface(226) = {225};
Line Loop(227) = {168, 73, -189, -191};
Ruled Surface(228) = {227};
Line Loop(229) = {166, 73, 190, -193};
Ruled Surface(230) = {229};
Line Loop(231) = {201, 210, -14, -212};
Ruled Surface(232) = {231};
Line Loop(233) = {15, 213, -202, 210};
Ruled Surface(234) = {233};
Line Loop(235) = {213, -200, 192, -171};
Ruled Surface(236) = {235};
Line Loop(237) = {172, -74, 188, 192};
Ruled Surface(238) = {237};
Line Loop(239) = {170, -74, -187, 194};
Ruled Surface(240) = {239};
Line Loop(241) = {187, -76, 28, -119};
Ruled Surface(242) = {241};
Line Loop(243) = {199, 212, 169, -194};
Ruled Surface(244) = {243};
Line Loop(245) = {119, 25, 75, 190};
Ruled Surface(246) = {245};
Line Loop(247) = {76, 188, 83, 27};
Ruled Surface(248) = {247};
Line Loop(249) = {26, -83, 189, -75};
Ruled Surface(250) = {249};
Line Loop(251) = {10, 11, 12, 9};
Line Loop(252) = {169, 170, -172, -171, -15, -14};
Plane Surface(253) = {251, 252};
Line Loop(254) = {41, 11, -43, -3};
Ruled Surface(255) = {254};
Line Loop(256) = {43, 12, -44, -4};
Ruled Surface(257) = {256};
Line Loop(258) = {44, 9, 42, -1};
Ruled Surface(259) = {258};
Line Loop(260) = {42, 2, 41, -10};
Ruled Surface(261) = {260};
Line Loop(262) = {3, 4, 1, 2};
Plane Surface(263) = {262};
Line Loop(264) = {38, 39, 40, 37};
Plane Surface(265) = {264};
Line Loop(266) = {37, -45, -33, 47};
Ruled Surface(267) = {266};
Line Loop(268) = {47, -40, 48, 36};
Ruled Surface(269) = {268};
Line Loop(270) = {48, -35, 46, 39};
Ruled Surface(271) = {270};
Line Loop(272) = {46, -38, -45, 34};
Ruled Surface(273) = {272};
Line Loop(274) = {33, 34, 35, 36};
Line Loop(275) = {181, 182, 29, 30, 179, 180};
Plane Surface(276) = {274, 275};
Line Loop(277) = {18, 207, 205, -217};
Ruled Surface(278) = {277};
Line Loop(279) = {207, -204, -215, -19};
Ruled Surface(280) = {279};
Line Loop(281) = {167, 197, 203, -215};
Ruled Surface(282) = {281};
Line Loop(283) = {217, 206, -195, -165};
Ruled Surface(284) = {283};
Line Loop(285) = {195, 183, -77, -166};
Ruled Surface(286) = {285};
Line Loop(287) = {197, -184, -77, -168};
Ruled Surface(288) = {287};
Line Loop(289) = {204, 208, 30, -216};
Ruled Surface(290) = {289};
Line Loop(291) = {203, 216, 179, -198};
Ruled Surface(292) = {291};
Line Loop(293) = {208, -29, -218, -205};
Ruled Surface(294) = {293};
Line Loop(295) = {198, 180, -78, -185};
Line Loop(296) = {218, -182, -196, -206};
Ruled Surface(297) = {296};
Line Loop(298) = {25, 113, -183, -85};
Ruled Surface(299) = {298};
Line Loop(300) = {85, -186, -114, 28};
Ruled Surface(301) = {300};
Line Loop(302) = {84, 185, -114, -27};
Ruled Surface(303) = {302};
Line Loop(304) = {26, 84, -184, -113};
Ruled Surface(305) = {304};
Ruled Surface(306) = {295};
Line Loop(307) = {78, 181, -196, -186};
Ruled Surface(308) = {307};
//Volume
Surface Loop(309) = {265, 273, 271, 269, 267, 276, 306, 292, 282, 226, 224, 280, 278, 222, 232, 234, 253, 261, 259, 257, 255, 263, 236, 238, 240, 242, 248, 250, 305, 303, 301, 299, 246, 230, 286, 284, 297, 294, 290, 308, 220, 244, 288, 228};
Volume(310) = {309};
//Physical Surface
Physical Surface("wall") = {273, 267, 269, 271, 276, 306, 292, 290, 282, 294, 288, 305, 308, 301, 286, 284, 280, 278, 297, 299, 303, 226, 222, 224, 228, 230, 246, 250, 220, 242, 248, 232, 234, 253, 236, 261, 238, 240, 259, 257, 255, 244};
Physical Surface("inlet") = {263};
Physical Surface("outlet") = {265};
//hysical volume
Physical Volume("volume") = {310};
