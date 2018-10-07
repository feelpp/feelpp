lc = 0.025;
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
Rotate {{0, 1, 0}, {0, 0, 0}, 2*Pi/3} {
  Duplicata { Line{287, 281, 278, 280, 279, 301, 297, 296, 305, 284, 293, 292, 291, 294}; }
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
Line Loop(396) = {353, -361, 352, -86, 339, -347, 338, -41, 281, -284, 287, -131};
Plane Surface(397) = {325, 396};
Line Loop(398) = {78, 354, 356, 355, 123, 278, 279, 280, 33, 340, 342, 341};
Line Loop(399) = {390, 391, 392};
Plane Surface(400) = {398, 399};
Line Loop(401) = {395, -391, -394, 318};
Ruled Surface(402) = {401};
Line Loop(403) = {394, -390, -393, 317};
Ruled Surface(404) = {403};
Line Loop(405) = {392, -393, -316, 395};
Ruled Surface(406) = {405};
Surface Loop(407) = {83, 371, 400, 369, 383, 367, 128, 132, 397, 331, 335, 333, 337, 404, 402, 406, 377, 310, 311, 302, 375, 38, 42, 373, 387, 381, 389, 87, 385, 379};
Volume(408) = {407};
Physical Surface(409) = {311, 381, 379};
Physical Surface(410) = {367, 383, 369, 385, 331, 333, 335, 397, 38, 42, 377, 310, 375, 302, 400, 132, 128, 83, 87, 371, 373, 389, 387};
Physical Surface(411) = {406, 404, 402, 337};
Physical Volume(412) = {408};
