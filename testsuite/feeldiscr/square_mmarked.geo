h=0.3;
Point(1) = {0, 1, 0, h};
Point(2) = {0, 2, 0, h};
Point(3) = {0, 3, 0, h};
Point(4) = {0, 4, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Extrude {1, 0, 0} {
  Line{1, 2, 3};
}
Extrude {1, 0, 0} {
  Line{12, 8, 4};
}
Extrude {1, 0, 0} {
  Line{16, 20, 24};
}
Extrude {0, 0, 1} {
  Surface{15, 19, 31, 35, 23, 11, 7, 27, 39};
}
Extrude {0, 0, 1} {
  Surface{61, 83, 105, 127, 149, 171, 193, 215, 237};
}
Extrude {0, 0, 1} {
  Surface{259, 281, 303, 369, 347, 325, 391, 413, 435};
}
Physical Volume("Omega") = {27, 24, 21, 18, 13, 12, 20, 3, 23, 4, 26, 9, 17, 14, 11, 19, 22, 25, 8, 16, 5, 2, 10, 15, 1, 6, 7};
Physical Point("P") = {110};
Physical Line("L") = {66};
Physical Surface("S") = {408};
Physical Surface("Wall1") = {48, 246, 444, 158, 356, 510, 180, 378, 576};
Physical Surface("Wall2") = {457, 479, 501, 523, 545, 567, 589, 611, 633};
Physical Surface("Wall3") = {496, 298, 100, 562, 320, 122, 628, 430, 232};
Physical Surface("Wall4") = {31, 19, 15, 35, 23, 11, 39, 27, 7};
Physical Surface("Wall5") = {96, 294, 492, 74, 272, 470, 448, 250, 52};
Physical Surface("Wall6") = {236, 214, 192, 434, 412, 390, 632, 610, 588};
