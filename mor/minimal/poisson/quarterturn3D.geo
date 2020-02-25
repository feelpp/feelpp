h = 0.1;

Point(1) = {1,0,0,h};
Point(2) = {2,0,0,h};

Line(1) = {1,2};
out[] = Extrude{{0,0,1},{0,0,0},Pi/2} {Line{1}; };
out1[] = Extrude{0,0,2} {Surface{out[1]}; };

Physical Volume("omega") = {out1[1]};
Physical Surface("bottom") = {out[1]};
Physical Surface("top") = {out1[0]};
Physical Surface("V1") = {out1[2]};
Physical Surface("V0") = {out1[4]};
Physical Surface("R") = {out1[3],out1[5]};
