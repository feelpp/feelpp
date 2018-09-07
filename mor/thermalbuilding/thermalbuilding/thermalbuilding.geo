h=0.1;

W=10;
w=8;
ep1=0.1;

Point(1) = {0,0,0.0,h};
Point(2) = {W,0,0.0,h};
Point(3) = {W,w,0.0,h};
Point(4) = {0,w,0.0,h};

Point(5) = {4,0,0.0,h};
Point(6) = {5,0,0.0,h};
Point(7) = {4-ep1,0,0.0,h};
Point(8) = {5+ep1,0,0.0,h};

Point(9) = {0,3.5-ep1,0.0,h};
Point(10) = {0,3.5,0.0,h};
Point(11) = {0,4.5,0.0,h};
Point(12) = {0,4.5+ep1,0.0,h};

Point(13) = {4-ep1,3.5-ep1,0.0,h};
Point(14) = {4,3.5,0.0,h};

Point(15) = {2,3.5,0.0,h};
Point(16) = {2,4.5,0.0,h};

Point(17) = {3,4.5,0.0,h};
Point(18) = {3-ep1,4.5+ep1,0.0,h};

Point(19) = {3-ep1,w,0.0,h};
Point(20) = {3,w,0.0,h};

Point(21) = {3,4.9,0.0,h};
Point(22) = {3-ep1,4.9,0.0,h};
Point(23) = {3,5.7,0.0,h};
Point(24) = {3-ep1,5.7,0.0,h};
Point(25) = {3,5.7+ep1,0.0,h};

Point(26) = {6,5.7,0.0,h};
Point(27) = {6-ep1,5.7+ep1,0.0,h};
Point(28) = {6-ep1,5.7,0.0,h};

Point(29) = {6,w,0.0,h};
Point(30) = {6-ep1,w,0.0,h};


Point(31) = {6,4.9,0.0,h};
Point(32) = {6-ep1,4.9,0.0,h};

//Point(33) = {6-ep1,4.6-ep1,0.0,h};
Point(34) = {6,4.6,0.0,h};
Point(35) = {6-ep1,4.6,0.0,h};

Point(36) = {5,4.6,0.0,h};
Point(37) = {5.1,4.6-ep1,0.0,h};

Point(38) = {5,3,0.0,h};
Point(39) = {5.1,3,0.0,h};

Point(40) = {5,1.5,0.0,h};
Point(41) = {5.1,1.5,0.0,h};

Point(42) = {W,4.6,0.0,h};
Point(43) = {W,4.6-ep1,0.0,h};

Point(44) = {4,3,0.0,h};
Point(45) = {4-ep1,3,0.0,h};

Point(46) = {4,1.5,0.0,h};
Point(47) = {4-ep1,1.5,0.0,h};

Point(48) = {2-ep1,3.5,0.0,h};
Point(49) = {2-ep1,4.5,0.0,h};


Point(50) = {3.5,5.7,0.0,h};
Point(51) = {3.5,5.7+ep1,0.0,h};
Point(52) = {4.3,5.7,0.0,h};
Point(53) = {4.3,5.7+ep1,0.0,h};

heaterHeightFoot=0.2;
Point(54) = {1.5,ep1,heaterHeightFoot,h};
Point(55) = {2.3,ep1,heaterHeightFoot,h};
Point(56) = {1.5,2*ep1,heaterHeightFoot,h};
Point(57) = {2.3,2*ep1,heaterHeightFoot,h};

Point(58) = {7,ep1,heaterHeightFoot,h};
Point(59) = {8,ep1,heaterHeightFoot,h};
Point(60) = {7,2*ep1,heaterHeightFoot,h};
Point(61) = {8,2*ep1,heaterHeightFoot,h};

Point(62) = {ep1,6,heaterHeightFoot,h};
Point(63) = {2*ep1,6,heaterHeightFoot,h};
Point(64) = {ep1,6.7,heaterHeightFoot,h};
Point(65) = {2*ep1,6.7,heaterHeightFoot,h};

Point(66) = {5.9-ep1,6,heaterHeightFoot,h};
Point(67) = {5.9-ep1,6.4,heaterHeightFoot,h};
Point(68) = {5.9-2*ep1,6,heaterHeightFoot,h};
Point(69) = {5.9-2*ep1,6.4,heaterHeightFoot,h};


Point(70) = {W-ep1,7,heaterHeightFoot,h};
Point(71) = {W-ep1,6.3,heaterHeightFoot,h};
Point(72) = {W-2*ep1,7,heaterHeightFoot,h};
Point(73) = {W-2*ep1,6.3,heaterHeightFoot,h};




Line(1) = {1, 7};
Line(2) = {7, 5};
Line(3) = {5, 6};
Line(4) = {6, 8};
Line(5) = {8, 2};
Line(6) = {2, 43};
Line(7) = {43, 42};
Line(8) = {42, 3};
Line(9) = {3, 29};
Line(10) = {29, 30};
Line(11) = {30, 20};
Line(12) = {20, 19};
Line(13) = {19, 4};
Line(14) = {4, 12};
Line(15) = {12, 11};
Line(16) = {11, 10};
Line(17) = {10, 9};
Line(18) = {9, 1};
Line(19) = {7, 47};
Line(20) = {5, 46};
Line(21) = {46, 47};
Line(22) = {45, 13};
Line(23) = {13, 9};
//Line(24) = {10, 15};
Line(25) = {15, 14};
Line(26) = {14, 44};
Line(27) = {44, 45};
Line(28) = {6, 40};
Line(29) = {40, 41};
Line(30) = {41, 8};
Line(31) = {38, 36};
Line(32) = {36, 35};
Line(33) = {35, 32};
Line(34) = {32, 31};
Line(35) = {31, 34};
Line(36) = {34, 42};
Line(37) = {43, 37};
Line(38) = {37, 39};
Line(39) = {39, 38};
Line(40) = {21, 22};
Line(41) = {22, 24};
Line(42) = {24, 23};
Line(43) = {23, 21};
Line(44) = {25, 20};
Line(45) = {18, 22};
Line(46) = {21, 17};
Line(47) = {17, 16};
//Line(48) = {16, 11};
Line(49) = {12, 18};
Line(50) = {24, 19};
Line(51) = {15, 48};
Line(52) = {48, 10};
Line(53) = {11, 49};
Line(54) = {49, 16};
Line(55) = {16, 15};
Line(56) = {48, 49};
Line(57) = {23, 50};
Line(58) = {50, 52};
Line(59) = {52, 28};
Line(60) = {25, 51};
Line(61) = {51, 53};
Line(62) = {53, 52};
Line(63) = {51, 50};
Line(64) = {53, 27};
Line(65) = {28, 26};
Line(66) = {26, 31};
Line(67) = {32, 28};
Line(68) = {26, 29};
Line(69) = {30, 27};
Line(70) = {54, 55};
Line(71) = {55, 57};
Line(72) = {57, 56};
Line(73) = {56, 54};
Line(74) = {58, 59};
Line(75) = {59, 61};
Line(76) = {61, 60};
Line(77) = {60, 58};
Line(78) = {71, 70};
Line(79) = {70, 72};
Line(80) = {72, 73};
Line(81) = {73, 71};
Line(82) = {69, 67};
Line(83) = {67, 66};
Line(84) = {66, 68};
Line(85) = {68, 69};
Line(86) = {65, 64};
Line(87) = {64, 62};
Line(88) = {62, 63};
Line(89) = {63, 65};


Line Loop(90) = {19, -21, -20, -2};
Plane Surface(91) = {90};
Line Loop(92) = {29, 30, -4, 28};
Plane Surface(93) = {92};
Line Loop(94) = {32, 33, 34, 35, 36, -7, 37, 38, 39, 31};
Plane Surface(95) = {94};
Line Loop(96) = {17, -23, -22, -27, -26, -25, 51, 52};
Plane Surface(97) = {96};
Line Loop(98) = {53, 54, -47, -46, 40, -45, -49, 15};
Plane Surface(99) = {98};
Line Loop(100) = {57, -63, -60, 44, 12, -50, 42};
Plane Surface(101) = {100};
Line Loop(102) = {64, -69, -10, -68, -65, -59, -62};
Plane Surface(103) = {102};

Line Loop(104) = {54, 55, 51, 56};
Plane Surface(105) = {104};
Line Loop(106) = {61, 62, -58, -63};
Plane Surface(107) = {106};
Line Loop(108) = {43, 40, 41, 42};
Plane Surface(109) = {108};
Line Loop(110) = {65, 66, -34, 67};
Plane Surface(111) = {110};

Line Loop(112) = {1, 19, -21, -20, 3, 28, 29, 30, 5, 6, 37, 38, 39, 31, 32, 33, 67, -59, -58, -57, 43, 46, 47, 55, 25, 26, 27, 22, 23, 18};
Line Loop(113) = {72, 73, 70, 71};
Line Loop(114) = {76, 77, 74, 75};
//Plane Surface(115) = {112, 113, 114};
Plane Surface(115) = {112};
Line Loop(116) = {53, -56, 52, -16};
Plane Surface(117) = {116};
Line Loop(118) = {13, 14, 49, 45, 41, 50};
Line Loop(119) = {89, 86, 87, 88};
//Plane Surface(120) = {118, 119};
Plane Surface(120) = {118};
Line Loop(121) = {64, -69, 11, -44, 60, 61};
Line Loop(122) = {85, 82, 83, 84};
//Plane Surface(123) = {121, 122};
Plane Surface(123) = {121};
Line Loop(124) = {68, -9, -8, -36, -35, -66};
Line Loop(125) = {80, 81, 78, 79};
//Plane Surface(126) = {124, 125};
Plane Surface(126) = {124};

heightBuilding=2.4;

Extrude {0, 0, heightBuilding} {
  Line{13, 12, 11, 9, 10, 8, 7, 6, 5, 4, 3, 2, 1, 18, 17, 16, 15, 14, 50, 44, 41, 43, 46, 45, 49, 47, 53, 54, 56, 55, 52, 23, 25, 26, 22, 27, 39, 31, 38, 32, 37, 33, 35, 36, 66, 67, 65, 59, 64, 62, 61, 60, 42, 40, 57, 58, 68, 69, 19, 21, 20, 28, 29, 30};
}
Extrude {0, 0, 0.5} {
  Line{87, 86, 89, 88, 82, 85, 83, 84, 79, 80, 78, 81, 72, 73, 70, 71, 77, 76, 75, 74};
}
Line(463) = {104, 106};
Line Loop(464) = {247, 183, -251, -263, -267, -259, -255, -463};
Plane Surface(465) = {464};
Line Loop(466) = {243, -463, 239, 235};
Plane Surface(467) = {466};
Line Loop(468) = {191, 231, 235, -227, -215, 339, -219, -223};
Plane Surface(469) = {468};
Line Loop(470) = {207, 335, 211, 339};
Plane Surface(471) = {470};
Line(472) = {132, 130};
Line Loop(473) = {131, -199, 335, 343, 472, -331, 203};
Plane Surface(474) = {473};
Line Loop(475) = {327, 323, -347, 472};
Plane Surface(476) = {475};
Line Loop(477) = {319, -355, -143, -351, -311, -315, -323};
Plane Surface(478) = {477};
Line(479) = {121, 122};
Line Loop(480) = {311, 303, -479, 307};
Plane Surface(481) = {480};
Line Loop(482) = {291, 479, 295, 299, -151, 287, 279, 271, 275, 283};
Plane Surface(483) = {482};
Line Loop(484) = {379, -163, 371, 375};
Plane Surface(485) = {484};
Line Loop(486) = {363, -359, 171, 367};
Plane Surface(487) = {486};
Line Loop(488) = {435, 439, 443, 431};
Plane Surface(489) = {488};
Line Loop(490) = {451, 447, 459, 455};
Plane Surface(491) = {490};
Line Loop(492) = {387, 383, 395, 391};
Plane Surface(493) = {492};
Line Loop(494) = {411, 403, 399, 407};
Plane Surface(495) = {494};
Line Loop(496) = {415, 419, 427, 423};
Plane Surface(497) = {496};
Plane Surface(498) = {113};
Plane Surface(499) = {114};
Plane Surface(500) = {125};
Plane Surface(501) = {119};
Plane Surface(502) = {122};
Line Loop(503) = {296, -479, -293, 34};
Plane Surface(504) = {503};
Line Loop(505) = {240, 463, -245, 51};
Plane Surface(506) = {505};
Line Loop(507) = {63, 345, 472, -328};
Plane Surface(508) = {507};
Surface Loop(509) = {362, 366, 487, 174, 91, 370};
Volume(510) = {509};
Surface Loop(511) = {374, 378, 382, 485, 166, 93};
Volume(512) = {511};
Surface Loop(513) = {254, 465, 250, 186, 97, 266, 270, 262, 258, 506};
Volume(514) = {513};
Surface Loop(515) = {234, 469, 194, 99, 230, 218, 222, 226, 238, 342};
Volume(516) = {515};
Surface Loop(517) = {202, 474, 134, 101, 346, 334, 206, 338, 508};
Volume(518) = {517};
Surface Loop(519) = {322, 358, 478, 146, 103, 354, 318, 326, 314};
Volume(520) = {519};
Surface Loop(521) = {302, 483, 294, 286, 278, 274, 282, 290, 298, 154, 95, 504};
Volume(522) = {521};

Surface Loop(523) = {467, 246, 242, 105, 238, 506};
Volume(524) = {523};
Surface Loop(525) = {214, 471, 210, 342, 338, 109};
Volume(526) = {525};
Surface Loop(527) = {508, 350, 476, 330, 326, 107};
Volume(528) = {527};
Surface Loop(529) = {314, 306, 481, 310, 504, 111};
Volume(530) = {529};



Surface Loop(531) = {489, 438, 442, 446, 434, 498};
Line Loop(532) = {179, 175, 359, -363, -367, 167, 371, 375, 379, 159, 155, 287, 279, 271, 275, 283, 291, 307, -315, -347, -343, 211, 215, 227, 243, 255, 259, 267, 263, 251};
Plane Surface(533) = {532};
Line Loop(534) = {247, -187, 231, -239};
Plane Surface(535) = {534};
Line Loop(536) = {223, 219, 207, 199, 127, 195};
Plane Surface(537) = {536};
Line Loop(538) = {203, -135, 355, -319, -327, -331};
Plane Surface(539) = {538};
Line Loop(540) = {139, -351, 303, 295, 299, 147};
Plane Surface(541) = {540};
Surface Loop(542) = {115, 178, 533, 182, 170, 162, 158, 434, 498, 438, 442, 446, 489, 370, 366, 362, 262, 270, 266, 254, 258, 282, 274, 278, 286, 294, 290, 318, 350, 214, 230, 218, 246, 310, 374, 378, 382, 454, 499, 450, 462, 458, 491, 346};

Volume(543) = {542};
Surface Loop(544) = {535, 190, 117, 242, 234, 250};
Volume(545) = {544};
Surface Loop(546) = {537, 130, 120, 198, 202, 226, 222, 394, 501, 390, 386, 398, 493, 210};
Volume(547) = {546};
Surface Loop(548) = {123, 138, 539, 206, 334, 358, 322, 330};
//Volume(549) = {548};
Surface Loop(550) = {126, 142, 541, 150, 354, 306, 298, 302};
Surface Loop(551) = {422, 500, 430, 426, 418, 497};
Volume(552) = {550, 551};

Surface Loop(553) = {138, 123, 539, 358, 322, 206, 334, 406, 502, 402, 410, 414, 495, 330};
Volume(549) = {553};

Physical Volume("internal-walls") = {510,512,514,516,518,520,522 };
Physical Volume("internal-doors") = { 524,526,528,530 };
Physical Volume("air") = {543,545,547,549,552 };

Physical Surface("front-door") = {170};
Physical Surface("exterior-walls") = {182, 178, 174, 170, 166, 162, 158, 154, 150, 142, 146, 138, 134, 130, 198, 194, 190, 186};
Physical Surface("heater-kitchen") = {489, 442, 434, 498, 438, 446};
Physical Surface("heater-livingroom") = {491, 450, 454, 499, 462, 458};
Physical Surface("heater-bedroom1") = {497, 422, 430, 418, 500, 426};
Physical Surface("heater-bedroom2") = {394, 386, 398, 390, 501, 493};
Physical Surface("heater-bathroom") = {495, 502, 402, 406, 410, 414};
