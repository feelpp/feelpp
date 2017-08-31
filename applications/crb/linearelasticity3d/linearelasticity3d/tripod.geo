meshIndex=1;
//lc = 0.025/(meshIndex+1);
lc = 0.015/(meshIndex+1);
offset = 0.075;
ratio = 1.;
Point(1) = {offset, 0.0, 0, lc};
Point(2) = {offset+0.1, 0.0, 0, lc};
Point(3) = {offset+0.1, 0.1, 0, lc};
Line(2) = {1, 2};
Line(3) = {2, 3};
Point(5) = {offset+0.08, 0.1, 0, lc*ratio};
Point(6) = {offset+0.08, 0.12, 0, lc*ratio};
Circle(4) = {3, 5, 6};
Point(7) = {offset+0.015, 0.14, 0, lc*ratio};
Point(8) = {offset+0.035, 0.12, 0, lc*ratio};
Point(9) = {offset+0.035, 0.14, 0, lc*ratio};
Point(10) = {offset, 0.14, 0, lc*ratio};
Line(5) = {6, 8};
Line(6) = {7, 10};
Line(7) = {10, 1};
Circle(8) = {8, 9, 7};
Line Loop(9) = {8, 6, 7, 2, 3, 4, 5};
Plane Surface(10) = {9};
angle = Pi/12;
Rotate {{0, 1, 0}, {0, 0, 0}, angle} {
      Surface{10};
}
Extrude {{0, 1, 0}, {0, 0, 0}, 2*Pi/3-2*angle} {
      Surface{10};
}
Rotate {{0, 1, 0}, {0, 0, 0}, 2*angle} {
      Duplicata { Surface{47}; }
}
Extrude {{0, 1, 0}, {0, 0, 0}, 2*Pi/3-2*angle} {
      Surface{48};
}
Rotate {{0, 1, 0}, {0, 0, 0}, 2*angle} {
      Duplicata { Surface{92}; }
}
Extrude {{0, 1, 0}, {0, 0, 0}, 2*Pi/3-2*angle} {
      Surface{93};
}
length = 0.3;
Point(173) = {length, 0, (offset+0.035)*Sin[3/2*angle], lc*length*3};
Point(174) = {length, 0, -(offset+0.035)*Sin[3/2*angle], lc*length*3};
Line(278) = {168, 173};
Line(279) = {173, 174};
Line(280) = {174, 2};
Translate {0, 0.12, 0} {
      Duplicata { Point{174}; }
}
Translate {0, 0.12, 0} {
      Duplicata { Point{173}; }
}
Line(281) = {6, 178};
Line(282) = {178, 174};
Line(283) = {179, 173};
Line(284) = {179, 178};
Line(287) = {179, 177};
Line Loop(288) = {283, 279, -282, -284};
Plane Surface(289) = {288};
height = 0.1;
Extrude {{0, 0, 1}, {length, -height, 0}, -Pi/2} {
  Surface{289};
}
Delete {
  Volume{4};
}
Delete {
  Surface{289, 298, 306};
}
Delete {
  Line{283, 282};
}
Delete {
  Volume{3, 1, 2};
}
Delete {
  Surface{137, 10, 93, 92, 47, 48};
}
Delete {
  Surface{124, 34, 79, 46, 136, 91, 26, 71, 116};
}
Delete {
  Line{97, 60, 105, 2, 108, 5, 100, 63, 55, 18, 15, 52, 50, 13, 103, 6, 95, 58};
}
Delete {
  Surface{112, 75, 67, 22, 30, 120};
}
Delete {
  Line{57, 94, 96, 59, 49, 51, 14, 12, 8, 7, 102, 104};
}
Circle(312) = {112, 49, 136};
Circle(313) = {108, 49, 132};
Circle(314) = {106, 47, 130};
Circle(315) = {84, 49, 13};
Circle(316) = {88, 49, 17};
Circle(317) = {82, 47, 11};
Circle(318) = {160, 49, 10};
Circle(319) = {156, 49, 7};
Circle(320) = {154, 47, 8};
Delete {
  Line{110, 320, 20, 21, 25, 319, 111, 115, 318, 312, 313, 314, 70, 78, 66, 65, 317, 315, 316};
}
Delete {
  Point{154, 106, 11, 13, 17, 156, 160, 112, 108};
}
Circle(316) = {10, 49, 136};
Circle(317) = {136, 49, 88};
Circle(318) = {88, 49, 10};
Circle(319) = {132, 49, 84};
Circle(320) = {84, 49, 7};
Circle(321) = {7, 49, 132};
Circle(322) = {130, 47, 82};
Circle(323) = {82, 47, 8};
Circle(324) = {8, 47, 130};
Line Loop(325) = {323, 324, 322};
Line Loop(326) = {320, 321, 319};
Circle(327) = {84, 83, 82};
Circle(328) = {132, 131, 130};
Circle(329) = {7, 9, 8};
Line Loop(330) = {323, -329, -320, 327};
Ruled Surface(331) = {330};
Line Loop(332) = {322, -327, -319, 328};
Ruled Surface(333) = {332};
Line Loop(334) = {328, -324, -329, 321};
Ruled Surface(335) = {334};
Line Loop(336) = {318, 316, 317};
Plane Surface(337) = {326, 336};
/*
Rotate {{0, 1, 0}, {0, 0, 0}, 2*Pi/3} {
  Duplicata { Line{287, 281, 278, 280, 279, 301, 297, 296, 305, 284, 293, 292, 291, 294,396,397,398}; }
}
Rotate {{0, 1, 0}, {0, 0, 0}, -2*Pi/3} {
  Duplicata { Line{287, 281, 278, 280, 279, 301, 297, 296, 305, 284, 293, 292, 291, 294}; }
}
Line Loop(366) = {355, 98, 99, 353, 360, -362, -357};
Plane Surface(367) = {366};
Line Loop(368) = {354, 358, -364, -359, 352, -62, -61};
Plane Surface(369) = {368};
Line Loop(370) = {341, 53, 54, 339, 346, -348, -343};
Plane Surface(371) = {370};
Line Loop(372) = {345, 350, -344, -340, 16, 17, -338};
Plane Surface(373) = {372};
Line Loop(374) = {280, 3, 4, 281, 305, -293, -301};
Plane Surface(375) = {374};
Line Loop(376) = {296, 291, -297, -278, 106, 107, -287};
Plane Surface(377) = {376};
Line Loop(378) = {364, 363, 362, 365};
Plane Surface(379) = {378};
Line Loop(380) = {348, 351, 350, 349};
Plane Surface(381) = {380};
Line Loop(382) = {357, -363, -358, 356};
Ruled Surface(383) = {382};
Line Loop(384) = {360, 365, -359, 361};
Ruled Surface(385) = {384};
Line Loop(386) = {346, 351, -345, 347};
Ruled Surface(387) = {386};
Line Loop(388) = {343, -349, -344, 342};
Ruled Surface(389) = {388};
*/

Delete {
  Line{74, 29, 119};
}
Delete {
  Point{116, 21, 164, 12, 155, 107};
}


Circle(390) = {140, 65, 92};
Circle(391) = {92, 65, 1};
Circle(392) = {1, 65, 140};
Line(393) = {136, 140};
Line(394) = {88, 92};
Line(395) = {10, 1};

/*
Line Loop(396) = {353, -361, 352, -86, 339, -347, 338, -41, 281, -284, 287, -131};
Plane Surface(397) = {325, 396};
Line Loop(398) = {78, 354, 356, 355, 123, 278, 279, 280, 33, 340, 342, 341};
Line Loop(399) = {390, 391, 392};
Plane Surface(400) = {398, 399};
*/

/*
Line Loop(469) = {395, -391, -394, 318};
Ruled Surface(470) = {469};
Line Loop(471) = {394, -390, -393, 317};
Ruled Surface(472) = {471};
Line Loop(473) = {392, -393, -316, 395};
Ruled Surface(474) = {473};
*/
/*Surface Loop(407) = {83, 371, 400, 369, 383, 367, 128, 132, 397, 331, 335, 333, 337, 404, 402, 406, 377, 310, 311, 302, 375, 38, 42, 373, 387, 381, 389, 87, 385, 379};
Volume(408) = {407};
*/

/*
Physical Surface("AppliedDisp") = {337};
//Physical Surface("FixedDisp") = {311, 381, 379};
Physical Surface("base1") = {311};
Physical Surface("base2") = {381};
Physical Surface("base3") = {379};

Physical Surface("flux") = {402,404,406};

Physical Volume("Material") = {408};
*/

Circle(396) = {6, 47, 177};
Circle(397) = {3, 81, 172};
Circle(398) = {2, 65, 168};
Rotate {{0, 1, 0}, {0, 0, 0}, 2*Pi/3} {
  Duplicata { Line{287, 281, 278, 280, 279, 301, 297, 296, 305, 284, 293, 292, 291, 294, 396, 397, 398}; }
}
Rotate {{0, 1, 0}, {0, 0, 0}, -2*Pi/3} {
  Duplicata { Line{287, 281, 278, 280, 279, 301, 297, 296, 305, 284, 293, 292, 291, 294, 396, 397, 398}; }
}

/*
Line Loop(399) = {396, -107, -397, 4};
Ruled Surface(400) = {399};
Line Loop(401) = {106, -397, -3, 398};
Ruled Surface(402) = {401};
Line Loop(403) = {396, -287, 284, -281};
Plane Surface(404) = {403};
Line Loop(405) = {278, 279, 280, 398};
Plane Surface(406) = {405};
Line Loop(407) = {107, -287, 296, 291, -297, -278, 106};
Plane Surface(408) = {407};
Line Loop(409) = {3, 4, 281, 305, -293, -301, 280};
Plane Surface(410) = {409};

Line Loop(433) = {413, -17, -414, 54};
Ruled Surface(434) = {433};
Line Loop(435) = {414, -16, -415, 53};
Ruled Surface(436) = {435};
Line Loop(437) = {400, -408, 399, -413};
Plane Surface(438) = {437};
Line Loop(439) = {401, 403, 402, 415};
Plane Surface(440) = {439};
//Line Loop(441) = {400, 407, -409, -404, 402, 53, 54};
//Plane Surface(442) = {441};
Line Loop(443) = {399, -17, -16, 401, 405, -411, -406};
Plane Surface(444) = {443};
Line Loop(445) = {409, 412, 411, 410};
Plane Surface(446) = {445};
Line Loop(447) = {408, 407, 412, -406};
Ruled Surface(448) = {447};
Line Loop(449) = {430, -62, -431, 99};
Ruled Surface(450) = {449};
Line Loop(451) = {431, -61, -432, 98};
Ruled Surface(452) = {451};
Line Loop(453) = {423, -429, -424, -425};
Ruled Surface(454) = {453};
Line Loop(455) = {416, -430, 417, -425};
Plane Surface(456) = {455};
Line Loop(457) = {62, -416, 423, 428, -422, -418, 61};
Plane Surface(458) = {457};
Line Loop(459) = {418, 420, 419, 432};
Plane Surface(460) = {459};
Line Loop(461) = {427, 426, 429, 428};
Plane Surface(462) = {461};
Line Loop(463) = {417, 424, -426, -421, 419, 98, 99};
Plane Surface(464) = {463};
Line Loop(465) = {420, 421, -427, -422};
Ruled Surface(466) = {465};
Line Loop(467) = {410, -404, -403, 405};
Ruled Surface(468) = {467};


Line Loop(469) = {395, -391, -394, 318};
Ruled Surface(470) = {469};
Line Loop(471) = {394, -390, -393, 317};
Ruled Surface(472) = {471};
Line Loop(473) = {392, -393, -316, 395};
Ruled Surface(474) = {473};
Line Loop(475) = {86, -430, 131, -396, 41, -413};
Plane Surface(476) = {325, 475};
Line Loop(477) = {390, 391, 392};
Line Loop(478) = {123, -398, 33, -415, 78, -432};
Plane Surface(479) = {477, 478};
Surface Loop(480) = {476, 331, 335, 333, 337, 472, 470, 474, 479, 83, 87, 38, 42, 128, 132, 434, 436, 450, 452, 404, 408, 310, 311, 302, 406, 410};
Surface Loop(481) = {42, 476, 331, 335, 333, 337, 472, 470, 474, 479, 83, 87, 38, 128, 132, 402, 400, 452, 450, 436, 434};
Volume(482) = {481};
Surface Loop(483) = {444, 438, 442, 448, 446, 468, 440, 436, 434};
Volume(484) = {483};
Surface Loop(485) = {460, 458, 456, 464, 454, 462, 466, 452, 450};
Volume(486) = {485};
Surface Loop(487) = {410, 404, 408, 310, 311, 302, 406, 402, 400};
Volume(488) = {487};


Physical Surface("base1") = {446};
Physical Surface("base2") = {462};
Physical Surface("base3") = {311};
Physical Surface("cylinder") = {470, 474, 472};
Physical Surface("support-top") = {337};

Physical Volume("mat-central") = {482};
Physical Volume("mat-base1") = {484};
Physical Volume("mat-base2") = {486};
Physical Volume("mat-base3") = {488};
*/
//+
//Line Loop(433) = {400, 407, -409, -404, 402, 53, 54};
//+
//Plane Surface(434) = {433};
//+
Line(435) = {173, 179};
//+
Line Loop(436) = {435, 296, 291, -297};
//+
Plane Surface(437) = {436};
//+
//+
Line(438) = {172, 173};
//+
Line Loop(439) = {107, -287, -435, -438};
//+
Ruled Surface(440) = {439};
//+
Line Loop(441) = {438, -278, 106};
//+
Plane Surface(442) = {441};
//+
Line(443) = {174, 178};
//+
Line(444) = {174, 3};
//+
Line(445) = {221, 212};
//+
Line(446) = {221, 29};
//+
Line(447) = {224, 217};
//+
Line(448) = {100, 224};
//+
Line(449) = {259, 250};
//+
Line(450) = {124, 259};
//+
Line(451) = {262, 255};
//+
Line(452) = {148, 262};
//+
Line Loop(453) = {293, -305, -443, 301};
//+
Plane Surface(454) = {453};
//+
Line Loop(455) = {280, 3, -444};
//+
Plane Surface(456) = {455};
//+
Line Loop(457) = {287, -396, 281, -284};
//+
Plane Surface(458) = {457};
//+
Line Loop(459) = {4, 281, -443, 444};
//+
Ruled Surface(460) = {459};
//+
Line Loop(461) = {396, -107, -397, 4};
//+
Ruled Surface(462) = {461};
//+
Line Loop(463) = {398, 106, -397, -3};
//+
Ruled Surface(464) = {463};
//+
Line Loop(465) = {280, 398, 278, 279};
//+
Plane Surface(466) = {465};
//+
Line Loop(467) = {16, -446, -401};
//+
Plane Surface(468) = {467};
//+
Line Loop(469) = {445, 406, 411, -405};
//+
Plane Surface(470) = {469};
//+
Line Loop(471) = {410, 409, 412, 411};
//+
Plane Surface(472) = {471};
//+
Line Loop(473) = {401, 403, 402, 415};
//+
Plane Surface(474) = {473};
//+
Line Loop(475) = {413, -399, 408, -400};
//+
Plane Surface(476) = {475};
//+
Line Loop(477) = {402, 53, 448};
//+
Plane Surface(478) = {477};
//+
Line Loop(479) = {404, 409, -407, -447};
//+
Plane Surface(480) = {479};
//+
Line Loop(481) = {406, -412, -407, -408};
//+
Ruled Surface(482) = {481};
//+
Line Loop(483) = {404, -410, -405, 403};
//+
Ruled Surface(484) = {483};
//+
Line Loop(485) = {415, 16, -414, -53};
//+
Ruled Surface(486) = {485};
//+
Line Loop(487) = {414, 17, -413, -54};
//+
Ruled Surface(488) = {487};
//+
Line Loop(489) = {448, 447, -400, -54};
//+
Ruled Surface(490) = {489};
//+
Line Loop(491) = {446, 17, -399, -445};
//+
Ruled Surface(492) = {491};
//+
Line Loop(493) = {427, 426, 429, 428};
//+
Plane Surface(494) = {493};
//+
Line Loop(495) = {422, -428, -423, -449};
//+
Plane Surface(496) = {495};
//+
Line Loop(497) = {426, -424, -451, 421};
//+
Plane Surface(498) = {497};
//+
Line Loop(499) = {61, 450, -418};
//+
Plane Surface(500) = {499};
//+
Line Loop(501) = {452, 419, 98};
//+
Plane Surface(502) = {501};
//+
Line Loop(503) = {416, -430, 417, -425};
//+
Plane Surface(504) = {503};
//+
Line Loop(505) = {432, 418, 420, 419};
//+
Plane Surface(506) = {505};
//+
Line Loop(507) = {423, -429, -424, -425};
//+
Ruled Surface(508) = {507};
//+
Line Loop(509) = {430, -62, -431, 99};
//+
Ruled Surface(510) = {509};
//+
Line Loop(511) = {431, -61, -432, 98};
//+
Ruled Surface(512) = {511};
//+
Line Loop(513) = {62, -416, -449, -450};
//+
Ruled Surface(514) = {513};
//+
Line Loop(515) = {99, 417, -451, -452};
//+
Ruled Surface(516) = {515};
//+
Line Loop(517) = {420, 421, -427, -422};
//+
Ruled Surface(518) = {517};
//+
Line Loop(519) = {394, -390, -393, 317};
//+
Ruled Surface(520) = {519};
//+
Line Loop(521) = {392, -393, -316, 395};
//+
Ruled Surface(522) = {521};
//+
Line Loop(523) = {395, -391, -394, 318};
//+
Ruled Surface(524) = {523};
//+
//+
Line Loop(525) = {86, -430, 131, -396, 41, -413};
//+
Plane Surface(526) = {325, 525};
//+
Line Loop(527) = {33, -415, 78, -432, 123, -398};
//+
Line Loop(528) = {392, 390, 391};
//+
Plane Surface(529) = {527, 528};
//+
Surface Loop(530) = {522, 529, 38, 42, 526, 331, 335, 333, 337, 520, 524, 132, 128, 87, 83, 512, 510, 488, 486, 462, 464};
//+
Volume(531) = {530}; //mat central
//+
Surface Loop(532) = {458, 440, 437, 310, 311, 302, 466, 456, 460, 454, 442, 464, 462};
//+
Volume(533) = {532};
//+
Surface Loop(534) = {468, 492, 476, 482, 470, 472, 484, 480, 490, 478, 474, 488, 486};
//+
Volume(535) = {534};
//+
Surface Loop(536) = {508, 496, 518, 506, 500, 514, 504, 516, 498, 494, 502, 512, 510};
//+
Volume(537) = {536};


Physical Volume("mat-central") = {531};
Physical Volume("mat-base1") = {533};
Physical Volume("mat-base2") = {535};
Physical Volume("mat-base3") = {537};

Physical Surface("base1") = {311};
Physical Surface("base2") = {472};
Physical Surface("base3") = {494};
Physical Surface("cylinder") = {520,522,524};
Physical Surface("support-top") = {337};
